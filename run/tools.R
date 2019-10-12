library(dplyr)



############################## FUNCTION TO GET STRATIFIED TRAIN/TEST SPLIT
## INPUT: 
  # data - original dataframe 
  # target_col - column name of target
  # start_pos_col- column name of start position of window
  # chr_col - column name of chromosome
  # feature_cols - vector of features names 
  # train_ratio - ratio of data to place at train dataset
  # seed - random seed for reproducibility
## RETURNS: all_data - list of: 
  # x_train - train features dataset 
  # y_train - train target (factor)
  # x_test - test features dataset
  # y_test - test target (factor)
## Returns data split stratified by target, chromosome and position in chromosome (divided on 4 parts)
##############################
get_train_test_split <- function(data, target_col, start_pos_col, chr_col, feature_cols,
                                 train_ratio, seed){
  
  chrs <- c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )
  data[chr_col] <- lapply(data[chr_col], as.character)
  
  # get position bins
  all_pos_map <- data.frame()
  
  for (chr in chrs){
    starts <- data[data[chr_col] == chr, start_pos_col]
    starts <- starts[order(starts)]
    bin_name <- cut(starts, breaks = 4, labels = c("first_bin", "second_bin", "third_bin", "last_bin"),
                    include.lowest = TRUE, right = FALSE)
    pos_map <- data.frame(start = starts, pos_bin = bin_name, chr = chr)
    names(pos_map) <- c(start_pos_col, "pos_bin", chr_col)
    all_pos_map <- rbind(all_pos_map, pos_map)
  }
  all_pos_map[chr_col] <- lapply(all_pos_map[chr_col], as.character)
  data <- data %>%
    inner_join(all_pos_map, by = c(chr_col, start_pos_col))
  
  # create strata column by chromosome and quantile of position in chromosome
  data['strata'] <- paste0(data[chr_col], "_", data['pos_bin'])
  
  # split data on examples with positive and negative target
  rownames(data) <- NULL
  pos_label_data <- data[data[target_col] == 1, ]
  neg_label_data <- data[data[target_col] == 0, ]
  
  # split in train_ratio proportion stratified by strata column
  set.seed(seed)
  tr_neg <- createDataPartition(neg_label_data$strata, p = train_ratio, list = FALSE)
  data_train <- neg_label_data[tr_neg,]
  data_test  <- neg_label_data[-tr_neg,]
  
  set.seed(seed)
  tr_pos <- createDataPartition(pos_label_data$strata, p = train_ratio, list = FALSE)
  data_train <- rbind(data_train, pos_label_data[tr_pos, ])
  data_test <- rbind(data_test, pos_label_data[-tr_pos, ])
  
  x_train <- data_train[feature_cols]
  x_test <- data_test[feature_cols]
  
  y_train <- as.factor(data_train[, target_col])
  levels(y_train) <- c("X0", "X1")
  y_test <- as.factor(data_test[, target_col])
  levels(y_test) <- c("X0", "X1")
  
  all_data <- list()
  all_data[[1]] <- x_train
  all_data[[2]] <- y_train
  all_data[[3]] <- x_test
  all_data[[4]] <- y_test
  
  names(all_data) <- c("x_train", "y_train", "x_test", "y_test")
  
  return(all_data)
}


############################## FUNCTION TO LIMIT FEATURES OUTLIERS
## INPUT: 
  # x_train - train features dataset 
  # x_test - test features dataset
  # features - column names of features to transform
  # iqr_multiplicator - multiplicator for interquartile range for determination of borders
## RETURNS: all_data - list of: 
  # train - transformed train features dataset 
  # test - transformed test features dataset
## Returns dataset with only "features" features, transformed by limiting outliers (borders 
## are fitted on train and applied to train and test)
##############################
limit_outliers <- function(x_train, x_test, features, iqr_multiplicator = 3){
  
  # fit borders
  borders <- apply(x_train[features], MARGIN = 2, FUN = function(x){
    data.frame(
      lower_border = as.numeric(quantile(x, 0.25)) - iqr_multiplicator * IQR(x),
      upper_border = as.numeric(quantile(x, 0.75)) + iqr_multiplicator * IQR(x)
    )
  })
  
  x_train_tr <- x_train[features]
  x_test_tr <- x_test[features]
  
  # transform train and test
  for (col in features){
    
    x_train_tr[x_train_tr[col] < as.numeric(borders[[col]]['lower_border']), col] <- 
      as.numeric(borders[[col]]['lower_border'])
    x_train_tr[x_train_tr[col] > as.numeric(borders[[col]]['upper_border']), col] <- 
      as.numeric(borders[[col]]['upper_border'])
    
    x_test_tr[x_test_tr[col] < as.numeric(borders[[col]]['lower_border']), col] <- 
      as.numeric(borders[[col]]['lower_border'])
    x_test_tr[x_test_tr[col] > as.numeric(borders[[col]]['upper_border']), col] <- 
      as.numeric(borders[[col]]['upper_border'])
  }
  
  result <- list()
  result[['train']] <- x_train_tr
  result[['test']] <- x_test_tr
  
  return(result)
}


############################## FUNCTION TO GET QUALITY METRICS FOR TRAIN AND TEST
## INPUT: 
  # train_pred - dataframe with columns "X1"(predicted probability of class "1") and "target" 
  # with levels "X0" and "X1" (class "1")
  # test_pred - dataframe with columns "X1"(predicted probability of class "1") and "target" 
  # with levels "X0" and "X1" (class "1")
  # model - caret model for which varImp function is applicable
## RETURNS: results - list of: 
  # roc_auc - ROC AUC on train and test 
  # importance - feature importance from model
  # recall - recall data for different probability quantiles 
##############################
get_model_quality <- function(train_pred, test_pred, model){
  
  # ROC AUC
  tr_auc <- as.numeric(auc(cases = train_pred$X1[train_pred$target == "X1"], controls = train_pred$X1[train_pred$target == "X0"],
                direction = "<"))
  te_auc <- as.numeric(auc(cases = test_pred$X1[test_pred$target == "X1"], controls = test_pred$X1[test_pred$target == "X0"],
                direction = "<"))
  out_df <- data.frame(tr_auc = tr_auc, te_auc = te_auc)
  
  # importance
  imp <- varImp(model, scale = FALSE)$importance
  names(imp) <- "importance"
  imp$feature <- rownames(imp)  
  rownames(imp) <- NULL
  
  
  # recall
  test_pred <- test_pred[order(test_pred$X1, decreasing=TRUE), ]
  test_pred$target <- unlist(lapply(test_pred$target, as.character))
  test_pred$target <- gsub(pattern = "X", replacement = "", x = test_pred$target)
  test_pred$target <- as.numeric(test_pred$target)
  test_pos <- sum(test_pred$target)
  
  quant <- 1 - seq(0.5, 0.95, 0.05)
  ind_q <- round(nrow(test_pred) * quant)
  ratio_data <- as.data.frame(cbind(quant, sapply(ind_q, function(x) sum(test_pred$target[1:x]) / test_pos)))
  names(ratio_data) <- c("quantile", "recall")
  
  results <- list()
  results[['roc_auc']] <- out_df
  results[['importance']] <- imp
  results[['recall']] <- ratio_data
  
  return(results)
}
