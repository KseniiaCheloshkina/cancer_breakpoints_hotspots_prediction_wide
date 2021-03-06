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

n_cores <- 3
registerDoParallel(n_cores)

set.seed(7)

win_len <- 100000
win_len_upper <- win_len * 10

# load features data
data_path <- "data/datasets/" 
data <- read.csv(
  paste0(data_path, "dataset_", format(win_len, scientific = FALSE), ".csv")
)
hsp_cols <- grep(x = names(data), pattern = "hsp", value = TRUE)
data <- data %>%
  select(-hsp_cols)

# load target data
target_path <- "data/target/"
q_target <- read.csv(
  paste0(target_path, "q_bkpt_", format(win_len, scientific = FALSE), ".csv")
)
q_bkpt_cols <- grep(x = names(q_target), pattern = "q", value = TRUE)

hsp_bkpt_target <- read.csv(
  paste0(target_path, "hsp_bkpt_", format(win_len, scientific = FALSE), ".csv")
)
hsp_bkpt_cols <- grep(x = names(hsp_bkpt_target), pattern = "hsp", value = TRUE)

# check the same target columns and exclude
dupl_cols <- duplicated(t(hsp_bkpt_target))
dupl_cols <- names(hsp_bkpt_target)[dupl_cols]
length(dupl_cols)

dupl_cols <- duplicated(t(q_target))
dupl_cols <- names(q_target)[dupl_cols]
length(dupl_cols)


# exclude targets with small number of positive examples
excluded_levels <- list()
excluded_levels[['100000']] <- c(
  "hsp_99.95._blood_red", "hsp_99.95._bone_red", "hsp_99.95._brain_red", "hsp_99.95._breast_red", "hsp_99.95._liver_red",
  "hsp_99.95._ovary_red", "hsp_99.95._pancreatic_red", "hsp_99.95._prostate_red", "hsp_99.95._skin_red",  "hsp_99.95._uterus_red",
  "hsp_99.99._blood_red", "hsp_99.99._bone_red", "hsp_99.99._brain_red" ,"hsp_99.99._breast_red", "hsp_99.99._liver_red",
  "hsp_99.99._ovary_red", "hsp_99.99._pancreatic_red", "hsp_99.99._prostate_red", "hsp_99.99._skin_red", "hsp_99.99._uterus_red")
excluded_levels[['1000000']] <- c(
  "hsp_99.9._blood_red", "hsp_99.9._bone_red", "hsp_99.9._brain_red", "hsp_99.9._breast_red", "hsp_99.9._liver_red",
  "hsp_99.9._ovary_red", "hsp_99.9._pancreatic_red", "hsp_99.9._prostate_red", "hsp_99.9._skin_red", "hsp_99.9._uterus_red", 
  "hsp_99.95._blood_red", "hsp_99.95._bone_red", "hsp_99.95._brain_red", "hsp_99.95._breast_red", "hsp_99.95._liver_red",
  "hsp_99.95._ovary_red", "hsp_99.95._pancreatic_red", "hsp_99.95._prostate_red", "hsp_99.95._skin_red",  "hsp_99.95._uterus_red",
  "hsp_99.99._blood_red", "hsp_99.99._bone_red", "hsp_99.99._brain_red", "hsp_99.99._breast_red", "hsp_99.99._liver_red",
  "hsp_99.99._ovary_red", "hsp_99.99._pancreatic_red", "hsp_99.99._prostate_red", "hsp_99.99._skin_red", "hsp_99.99._uterus_red")
excl_lev_current <- excluded_levels[[as.character(format(win_len, scientific = F))]]
if (!is.null(excl_lev_current)){
  hsp_bkpt_cols <- setdiff(hsp_bkpt_cols, excl_lev_current)
}

all_target_cols <- c(hsp_bkpt_cols, q_bkpt_cols)


# join features and target
all_target_data <- q_target %>%
  inner_join(
    hsp_bkpt_target, by=c("chr", "from", "to")
  ) %>%
  select(c("chr", "from", "to", all_target_cols))

data <- data %>%
  inner_join(
    all_target_data, by=c("chr", "from", "to")
)



output_path <- "data/output/randomness/"

ss_cols <- c ("chr", "from", "to")
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

tissue_spec_feats <- setdiff(names(data), c(ss_cols, all_target_cols, conserved_features))

features_cols <- c(conserved_features, tissue_spec_feats)

all_data <- get_higher_level_features(data=data, features_cols = features_cols, 
                                           win_len_upper = win_len_upper,
                                           path_to_upper_data = paste0(data_path, "dataset_",
                                                                       format(win_len_upper, scientific = FALSE), ".csv"))

high_features <- setdiff(names(all_data), c(features_cols, ss_cols, all_target_cols))
features_cols <- c(features_cols, high_features)

all_tissue_spec_feats <- vector()
for (feat in features_cols){
  for (col in tissue_spec_feats){
    all_tissue_spec_feats <- c(all_tissue_spec_feats, grep(x = feat, pattern = col, value = TRUE))  
  }
}
all_conserved_feats <- setdiff(features_cols, all_tissue_spec_feats)

# train/test split
train_ratio <- 0.7
n_repeats <- 30

# for storing results
df_roc_auc_all <- data.frame()
df_recall_all <- data.frame()

# select quantiles for quality assessment
recall_quantiles <- c(0.001, 0.005, 0.002,  0.003, 0.015, 0.025, seq(0.01, 0.05, 0.01), seq(0.1, 0.9, 0.05))

# set progress bar
n_iterations <- length(all_target_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)

for (target_column in all_target_cols){

  cancer_type <- strsplit(target_column, "_")[[1]][3]
  agg_level <- strsplit(target_column, "_")[[1]][2]
    
  all_features_cols <- c(all_conserved_feats, grep(x = all_tissue_spec_feats,
                                                   pattern = cancer_type, value = TRUE))

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

  # if it is classification task for "breakpoints" vs "hotspots" then drop rows without any breakpoints
  if (nrow(unique(all_data[target_column])) == 3){
    new_all_data <- all_data[all_data[target_column] != 2, ]
  } else {
    new_all_data <- all_data
  }
    
  repeats_res <- foreach(i=seq(1, n_repeats)) %dopar% {

    # train/test split
    splitted_dataset <- get_train_test_split(data=new_all_data, target_col=target_column, 
                                             start_pos_col="from",
                                             chr_col="chr", feature_cols=all_features_cols,
                                             train_ratio=train_ratio, 
                                             seed=i)
    x_train <- splitted_dataset[["x_train"]]
    y_train <- splitted_dataset[["y_train"]]
    x_test <- splitted_dataset[["x_test"]] 
    y_test <- splitted_dataset[["y_test"]]  

    n_pos <- length(y_train[y_train == "X1"])
    n_neg <- length(y_train[y_train == "X0"])

    # fit models
    model <- rf_fit(x_train = x_train[all_features_cols], y_train = y_train, trCtrl = trCtrl,
                    n_pos = n_pos, n_neg = n_neg, mtryGrid = mtryGrid)

    # prediction
    train_pred <- predict(model, newdata = x_train[all_features_cols], type = "prob")
    test_pred <- predict(model, newdata = x_test[all_features_cols], type = "prob")
    train_pred$target <- y_train
    test_pred$target <- y_test

    # model quality
    model_qual <- get_model_quality(train_pred, test_pred, model, recall_quantiles)
    return(model_qual)
  }


  # save results
  for (split_iter in 1:length(repeats_res)){

    res_iter <- repeats_res[[split_iter]]

    # extract specific datasets
    df_roc_auc <- res_iter[['roc_auc']]
    df_recall <- res_iter[['recall']]

    # add general info
    df_roc_auc <- df_roc_auc %>%
      mutate(
        iter = split_iter,
        cancer_type = cancer_type,
        agg_level = agg_level,
        win_len = format(win_len, scientific = F)
      )

    df_recall <- df_recall %>%
      mutate(
        iter = split_iter,
        cancer_type = cancer_type,
        agg_level = agg_level,
        win_len = format(win_len, scientific = F)
      )

    df_roc_auc_all <- rbind(df_roc_auc_all, df_roc_auc)
    df_recall_all <- rbind(df_recall_all, df_recall)
    
    
    write.csv(df_roc_auc_all, file = paste0(output_path, "result_roc_auc_", format(win_len, scientific = F), ".csv"))
    write.csv(df_recall_all, file = paste0(output_path, "result_recall_", format(win_len, scientific = F), ".csv"))
  }  
  pb$tick()
}
