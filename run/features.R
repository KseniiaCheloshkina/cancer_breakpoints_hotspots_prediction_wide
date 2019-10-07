library(dplyr)

############################## FUNCTION TO CONVERT COVERAGE INTO BINARY STATE
## INPUT: data - dataframe for which each column will be transformed into binary (1/0) state
## RETURNS: data - dataframe: initial dataframe + new features
##############################
get_binary_features <- function(data, features_cols){
  new_features <- data.frame(apply(X = data[features_cols], MARGIN =  2, FUN = function(x) ifelse(x > 0, 1, 0)))
  new_names <- paste0("bin_", names(new_features))
  names(new_features) <- new_names  
  all_features <- cbind(data, new_features)
  return(all_features)
}



############################## FUNCTION TO GET HIGHER-LEVEL OF AGGREGATION FEATURES
## INPUT: 
    # data - original dataframe
    # win_len_upper - level of aggregation for which to take new (higher-level) features
    # path_to_upper_data - full path to datafrom higher level
    # features_cols - vector of all columns names which to take from higher-level data

## RETURNS: data - dataframe: initial dataframe + new features
##############################
get_higher_level_features <- function(data, win_len_upper, path_to_upper_data, features_cols){
  
  data_upper <- read.csv(path_to_upper_data)
  data_upper$to_window <- ceiling(data_upper$to / win_len_upper)
  
  data_upper <- data_upper %>%
    select(c("chr", "to_window", features_cols))
  
  nm_feats <- names(data_upper)[3:length(names(data_upper))]
  nm_feats <- paste0("upper_", nm_feats)
  names(data_upper)[3:length(names(data_upper))] <- nm_feats
  
  
  data$to_window <- ceiling(data$to / win_len_upper)
  new_features <- data %>%
    left_join(data_upper, by=c("chr", "to_window")) %>%
    select(-to_window)
  
  for (col in nm_feats){
    
    new_features <- new_features %>% 
      mutate(!!as.name(col) := if_else(is.na(!!as.name(col)), 0, !!as.name(col)))    
  }

  return(new_features)
}


############################## FUNCTION TO GET FLAGS OF WINDOW BEING LOCAL OR GLOBAL MAXIMUM
## INPUT: 
# data - original dataframe
# features_cols - vector of all columns for which to get flags
## RETURNS: data - dataframe: initial dataframe + new features
##############################
get_maximum_features <- function(data, features_cols){
  
  # global maximums
  q <- c(0.9, 0.95, 0.99)
  q_res <- apply(data[features_cols], MARGIN = 2, FUN = function(x) quantile(x, q))
  q_nm <- row.names(q_res)
  for (col in features_cols){
    for (nm in q_nm){
      new_col_name <- paste0("max_", nm, "_", col)
      data[new_col_name] <- ifelse(data[col] < q_res[nm, col], 0, 1)
    }
  }
  
  # local maximums
  data <- data %>%
    arrange(chr, from)
  
  for (col in features_cols){
    
    for (i in 1:10){
      data <- data %>% 
        mutate(
          !!paste0("lag_", i) := lag(!!as.name(col), n = i),
          !!paste0("lead_", i) := lead(!!as.name(col), n = i)
        )
    }
    
    # is element maximal for nearest 1, 5 or 10 neighbours (previous and following neighbours are considered)
    for (n_neighbours in c(1, 5, 10)){
      data$max_neighb <- apply(data[c(col, paste0("lag_", seq(1, n_neighbours)), paste0("lead_", seq(1, n_neighbours)))], 
                               MARGIN = 1, function(x) max(x, na.rm = TRUE))
      max_col_nm <- paste0("max_", n_neighbours, "_", col)
      data[max_col_nm] <- ifelse(data[col] == data['max_neighb'], 1, 0)  
    }
    bin_prefixes <- c("max_5_", "max_10_", "max_1_", "max_90%_", "max_95%_", "max_99%_")
    data[, paste0(bin_prefixes, col)] <- apply(
      data[, paste0(bin_prefixes, col)], MARGIN = 2,FUN = function(x) ifelse(is.na(x), 0, x))
    
    # max relative difference between element and 10 neighbours
    data$max_neighb <- apply(data[c(paste0("lag_", seq(1, 10)), paste0("lead_", seq(1, 10)))], 
                             MARGIN = 1, function(x) max(x, na.rm = TRUE))
    ratio_col_nm <- paste0("max_ratio_10_", col)
    data[ratio_col_nm] <- (data[col] - data['max_neighb']) / data[col] 
    
    lower_border <- quantile(data[is.finite(data[, ratio_col_nm]), ratio_col_nm], 0.25) - 
      3 * IQR(data[is.finite(data[, ratio_col_nm]), ratio_col_nm])
    data[is.na(data[ratio_col_nm]), ratio_col_nm] <- 0
    data[is.infinite(data[, ratio_col_nm]), ratio_col_nm] <- lower_border
    data[data[ratio_col_nm] < lower_border, ratio_col_nm] <- lower_border
    
    # drop temp columns
    data[c("max_neighb", paste0("lag_", seq(1, 10)), paste0("lead_", seq(1, 10)))] <- NULL
    
  }
  
  return(data)  
}


