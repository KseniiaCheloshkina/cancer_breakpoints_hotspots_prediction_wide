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

get_feature_df <- function(features_cols){

  # get groups of features
  sec_str <- c("A_Phased_Repeat","Direct_Repeat", "Inverted_Repeat", "Mirror_Repeat","Short_Tandem_Repeat",
               "G_quadruplex", "stemloops_16_50", "stemloops_6_15", "Z_DNA_Motif")
  all_sec_str <- vector()
  for (feat in features_cols){
    for (col in sec_str){
      all_sec_str <- c(all_sec_str, grep(x = feat, pattern = col, value = TRUE))  
    }
  }
  
  reg <- c("X3UTR", "X5UTR", "codingExons", "downstream", "introns", "promoters", "WholeGenes")
  all_reg <- vector()
  for (feat in features_cols){
    for (col in reg){
      all_reg <- c(all_reg, grep(x = feat, pattern = col, value = TRUE))  
    }
  }
  
  tad <- c("tad_boundaries_liver", "tad_boundaries_ovary", "tad_boundaries_pancreatic")
  all_tad <- vector()
  for (feat in features_cols){
    for (col in tad){
      all_tad <- c(all_tad, grep(x = feat, pattern = col, value = TRUE))  
    }
  }
  
  chromatin <-c("cancer_skin_DNase_seq", "cancer_brain_DNase_seq", "cancer_blood_DNase_seq",
                "cancer_prostate_DNase_seq", "cancer_ovary_DNase_seq",
                "cancer_liver_DNase_seq", "cancer_breast_DNase_seq", "cancer_uterus_DNase_seq",
                "cancer_bone_DNase_seq")
  all_chromatin <- vector()
  for (feat in features_cols){
    for (col in chromatin){
      all_chromatin <- c(all_chromatin, grep(x = feat, pattern = col, value = TRUE))  
    }
  }
  
  methyl <- c("cancer_brain_DNA_methylation", "cancer_breast_DNA_methylation",
              "cancer_liver_DNA_methylation",
              "cancer_skin_DNA_methylation", "cancer_uterus_DNA_methylation")
  all_methyl <- vector()
  for (feat in features_cols){
    for (col in methyl){
      all_methyl <- c(all_methyl, grep(x = feat, pattern = col, value = TRUE))  
    }
  }
  
  histones <- c( "cancer_liver_H3K27ac.human", "cancer_uterus_H3K27ac.human", "cancer_blood_H3K27ac.human",
                 "cancer_brain_H3K27ac.human","cancer_blood_H3K27me3.human","cancer_brain_H3K27me3.human",    
                 "cancer_blood_H3K36me3.human","cancer_brain_H3K36me3.human",    
                 "cancer_blood_H3K4me1.human","cancer_brain_H3K4me1.human",
                 "cancer_breast_H3K4me3.human","cancer_uterus_H3K4me3.human","cancer_liver_H3K4me3.human",
                 "cancer_brain_H3K4me3.human","cancer_blood_H3K4me3.human","cancer_skin_H3K4me3.human",
                 "cancer_brain_H3K9me3.human",
                 "cancer_blood_H3K9me3.human" )
  all_histones <- vector()
  for (feat in features_cols){
    for (col in histones){
      all_histones <- c(all_histones, grep(x = feat, pattern = col, value = TRUE))  
    }
  }
  
  tf <- c("cancer_liver_ATF3.human", "cancer_liver_CTCF.human",
          "cancer_liver_EGR1.human", "cancer_liver_FOXA1.human", "cancer_liver_FOXA2.human", 
          "cancer_liver_GABPA.human","cancer_liver_HNF4A.human", "cancer_liver_HNF4G.human",
          "cancer_liver_JUND.human","cancer_liver_MAX.human", "cancer_liver_NR2F2.human",
          "cancer_liver_REST.human", "cancer_liver_RXRA.human", "cancer_liver_SP1.human",
          "cancer_liver_YY1.human", "cancer_liver_ZBTB33.human")
  all_tf <- vector()
  for (feat in features_cols){
    for (col in tf){
      all_tf <- c(all_tf, grep(x = feat, pattern = col, value = TRUE))  
    }
  }
  feat_group_df <- rbind(
    data.frame(
      feature_group=c("tf"),
      color=c("#8F7700FF"),
      feature=all_tf
    ),
    data.frame(
      feature_group=c("histones"),
      color=c("#EFC000FF"),
      feature=all_histones
    ),
    data.frame(
      feature_group=c("methyl"),
      color=c("#868686FF"),
      feature=all_methyl
    ),
    data.frame(
      feature_group=c("chromatin"),
      color=c("#0073C2FF"),
      feature=all_chromatin
    ),
    data.frame(
      feature_group=c("tad"),
      color=c("#003С67FF"),
      feature=all_tad
    ),
    data.frame(
      feature_group=c("reg"),
      color=c("#7AA6DCFF"),
      feature=all_reg
    ),
    data.frame(
      feature_group=c("sec_str"),
      color=c("#A73030FF"),
      feature=all_sec_str
    )
  )
  return(feat_group_df)
}