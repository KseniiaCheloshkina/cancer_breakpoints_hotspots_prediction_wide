script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(reshape2)

library(progress)
library(doParallel)

library(caret)
library(randomForest)
library(pROC)

source("run/features.R")
source("run/tools.R")

n_cores <- 2
registerDoParallel(n_cores)

set.seed(7)


win_len <- 100000
win_len_upper <- win_len * 10

# load data
data_path <- "data/datasets/" 
data <- read.csv(
  paste0(data_path, "dataset_", format(win_len, scientific = FALSE), ".csv")
)
output_path <- "data/output/classifier_results_by_feature_group_predictions/" 

hsp_cols <- grep(x = names(data), pattern = "hsp", value = TRUE)

conserved_features <- c(
  "A_Phased_Repeat", "Direct_Repeat", 
  "Inverted_Repeat", "Mirror_Repeat", "Short_Tandem_Repeat", 
  "Z_DNA_Motif", "G_quadruplex", 
  "stemloops_16_50", "stemloops_6_15", 
  "X3UTR", "X5UTR", "codingExons", "downstream", "introns", "promoters", "WholeGenes", 
  "tad_boundaries_liver", "tad_boundaries_ovary", "tad_boundaries_pancreatic",
  
  "cancer_liver_ATF3.human", "cancer_liver_CTCF.human", "cancer_liver_EGR1.human",
  "cancer_liver_FOXA1.human", "cancer_liver_FOXA2.human", "cancer_liver_GABPA.human",
  "cancer_liver_HNF4A.human", "cancer_liver_HNF4G.human", "cancer_liver_JUND.human",
  "cancer_liver_MAX.human", "cancer_liver_NR2F2.human", "cancer_liver_REST.human",
  "cancer_liver_RXRA.human", "cancer_liver_SP1.human", "cancer_liver_YY1.human",
  "cancer_liver_ZBTB33.human"
)
tissue_spec_feats <- setdiff(names(data), c("chr","from","to", hsp_cols, conserved_features))

ss_cols <- c ("chr", "from", "to")
features_cols <- c(conserved_features, tissue_spec_feats)

all_data <- get_higher_level_features(data=data, features_cols = features_cols, 
                                      win_len_upper = win_len_upper,
                                      path_to_upper_data = paste0(data_path, "dataset_",
                                                                  format(win_len_upper, scientific = FALSE), ".csv"))

high_features <- setdiff(names(all_data), c(features_cols, ss_cols, hsp_cols))
features_cols <- c(features_cols, high_features)


all_tissue_spec_feats <- vector()
for (feat in features_cols){
  for (col in tissue_spec_feats){
    all_tissue_spec_feats <- c(all_tissue_spec_feats, grep(x = feat, pattern = col, value = TRUE))  
  }
}
all_conserved_feats <- setdiff(features_cols, all_tissue_spec_feats)

# get groups of features
df_features <- get_feature_df(features_cols)

# exclude targets with small number of positive examples
excluded_levels <- list()
excluded_levels[['100000']] <- c(
  "hsp_99.95._blood", "hsp_99.95._bone", "hsp_99.95._brain", "hsp_99.95._breast", "hsp_99.95._liver",
  "hsp_99.95._ovary", "hsp_99.95._pancreatic", "hsp_99.95._prostate", "hsp_99.95._skin",  "hsp_99.95._uterus",
  "hsp_99.99._blood", "hsp_99.99._bone", "hsp_99.99._brain", "hsp_99.99._breast", "hsp_99.99._liver",
  "hsp_99.99._ovary", "hsp_99.99._pancreatic", "hsp_99.99._prostate", "hsp_99.99._skin", "hsp_99.99._uterus")
excluded_levels[['1000000']] <- c(
  "hsp_99.9._blood", "hsp_99.9._bone", "hsp_99.9._brain", "hsp_99.9._breast", "hsp_99.9._liver",
  "hsp_99.9._ovary", "hsp_99.9._pancreatic", "hsp_99.9._prostate", "hsp_99.9._skin", "hsp_99.9._uterus", 
  "hsp_99.95._blood", "hsp_99.95._bone", "hsp_99.95._brain", "hsp_99.95._breast", "hsp_99.95._liver",
  "hsp_99.95._ovary", "hsp_99.95._pancreatic", "hsp_99.95._prostate", "hsp_99.95._skin",  "hsp_99.95._uterus",
  "hsp_99.99._blood", "hsp_99.99._bone", "hsp_99.99._brain", "hsp_99.99._breast", "hsp_99.99._liver",
  "hsp_99.99._ovary", "hsp_99.99._pancreatic", "hsp_99.99._prostate", "hsp_99.99._skin", "hsp_99.99._uterus")

excl_lev_current <- excluded_levels[[as.character(format(win_len, scientific = F))]]
if (!is.null(excl_lev_current)){
  hsp_cols <- setdiff(hsp_cols, excl_lev_current)
}


# train/test split
train_ratio <- 0.7
n_repeats <- 30


# for storing results
df_predictions <- data.frame()

# select quantiles for quality assessment
breakpoints_balance <- data.frame(t(apply(all_data[grep(x = hsp_cols, pattern = "all", value = T)], 2, table)))
breakpoints_balance$ratio <- breakpoints_balance$X1/ (breakpoints_balance$X0 + breakpoints_balance$X1)
recall_quantiles <- c(0.001, 0.005, 0.002,  0.003, 0.015, 0.025, seq(0.01, 0.05, 0.01), seq(0.1, 0.9, 0.05))

# predict only hotspots
hsp_cols <- hsp_cols[-grep(hsp_cols, pattern = "all")]
hsp_cols <- c("hsp_99.5._skin", "hsp_99.9._skin", "hsp_99._uterus", "hsp_99.5._uterus", "hsp_99.9._uterus")

# set progress bar
n_iterations <- length(hsp_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)


for (target_column in hsp_cols){
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  agg_level <- strsplit(target_column, "_")[[1]][2]
  
  all_features_cols <- c(all_conserved_feats, grep(x = all_tissue_spec_feats, 
                                                   pattern = cancer_type, value = TRUE))

  # define features list for iterative model evaluation
  feat_groups <- list()
  feat_groups[["sec_str"]] <- intersect(
    df_features[df_features['feature_group'] == 'sec_str', 'feature'] %>% as.character(), 
    all_features_cols
    )
  feat_groups[["reg"]] <- intersect(
    df_features[df_features['feature_group'] == 'reg', 'feature'] %>% as.character(), 
    all_features_cols
    )
  feat_groups[["tad"]] <- intersect(
    df_features[df_features['feature_group'] == 'tad', 'feature'] %>% as.character(), 
    all_features_cols
    )
  feat_groups[["tf"]] <- intersect(
    df_features[df_features['feature_group'] == 'tf', 'feature'] %>% as.character(), 
    all_features_cols
    )
  # cancer-specific
  int <- intersect(
    df_features[df_features['feature_group'] == 'histones', 'feature'] %>% as.character(), 
    all_features_cols
    )
  if (length(int) > 0){
    feat_groups[["histones"]] <- int
  }
  int <- intersect(
    df_features[df_features['feature_group'] == 'methyl', 'feature'] %>% as.character(),
    all_features_cols)
  if (length(int) > 0){
    feat_groups[["methyl"]] <- int
  }
  int <- intersect(
    df_features[df_features['feature_group'] == 'chromatin', 'feature'] %>% as.character(),
    all_features_cols
    )
  if (length(int) > 0){
    feat_groups[["chromatin"]] <- int
  }  

  # models parameters
  trCtrl <- trainControl(
    method="none",
    verboseIter = TRUE,
    classProbs=TRUE,
    summaryFunction = twoClassSummary,
    seeds = seq(1, n_repeats),
    allowParallel = TRUE
  )
  # set "seeds" parameters to make results reproducible
  mtryGrid <- expand.grid(mtry = 5)
  
  repeats_res <- foreach(i=seq(1, n_repeats)) %dopar% {
    
    # train/test split
    splitted_dataset <- get_train_test_split(data=all_data, target_col=target_column, 
                                             start_pos_col="from", chr_col="chr", 
                                             feature_cols=all_features_cols, train_ratio=train_ratio, 
                                             seed=i)
    x_train <- splitted_dataset[["x_train"]]
    y_train <- splitted_dataset[["y_train"]]
    x_test <- splitted_dataset[["x_test"]] 
    y_test <- splitted_dataset[["y_test"]]  

    feat_iter_res <- data.frame()
    # for each feature group
    for (j in 1:length(feat_groups)){
      
      features_gr <- names(feat_groups)[j]
      features_nm <- feat_groups[[features_gr]]
      
      n_pos <- length(y_train[y_train == "X1"])
      n_neg <- length(y_train[y_train == "X0"])
      
      # fit models
      model <- rf_fit(x_train = x_train[features_nm], y_train = y_train, trCtrl = trCtrl,
                      n_pos = n_pos, n_neg = n_neg, mtryGrid = mtryGrid)
      
      # prediction
      train_pred <- predict(model, newdata = x_train[features_nm], type = "prob")
      test_pred <- predict(model, newdata = x_test[features_nm], type = "prob")
      train_pred$target <- y_train
      test_pred$target <- y_test

      test_pred$X0 <- NULL
      test_pred$features_group <- features_gr
      test_pred$rn <- rownames(test_pred)
      
      feat_iter_res <- rbind(feat_iter_res, test_pred)
    }
    
    return(feat_iter_res)
  }
  
  
  # save results
  for (split_iter in 1:length(repeats_res)){
    res_iter <- repeats_res[[split_iter]]
    
    # add general info
    res_iter <- res_iter %>%
      mutate(
        repeats = split_iter,
        cancer_type = cancer_type,
        agg_level = agg_level,
        win_len = format(win_len, scientific = F)
      )
    
    df_predictions <- rbind(df_predictions, res_iter)  
  }
  
  write.csv(df_predictions, file = paste0(output_path, "result_pred_", format(win_len, scientific = F), ".csv"))
  pb$tick()
}



