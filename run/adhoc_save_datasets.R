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


############################################ HOTSPOTS CLASSIFICATION DATA AND BREAKPOINTS CLASSIFICATION DATA

# load data
data_path <- "data/datasets/classification/" 
data <- read.csv(
  paste0(data_path, "dataset_", format(win_len, scientific = FALSE), ".csv")
)

output_path <- "data/adhoc/datasets/"

# select only needed cols
hsp_cols <- grep(x = names(data), pattern = "hsp", value = TRUE)
ss_cols <- c ("chr", "from", "to")
other_cols <- setdiff(names(data), c(hsp_cols, ss_cols))
data[other_cols] <- NULL

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
# df_data <- data.frame()

# set progress bar
n_iterations <- length(hsp_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)

# target_column <- hsp_cols[1]
for (target_column in hsp_cols){
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  agg_level <- strsplit(target_column, "_")[[1]][2]
  
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
  
  repeats_res <- foreach(i=seq(1, n_repeats),.combine=rbind) %dopar% {
    
    # train/test split
    splitted_dataset <- get_train_test_split(data=data, target_col=target_column, 
                                             start_pos_col="from", chr_col="chr", 
                                             feature_cols=ss_cols, train_ratio=train_ratio, 
                                             seed=i)
    x_train <- splitted_dataset[["x_train"]]
    y_train <- splitted_dataset[["y_train"]]
    x_test <- splitted_dataset[["x_test"]] 
    y_test <- splitted_dataset[["y_test"]]
    
    # save into one dataset
    x_train['data_type'] <- 'train'
    x_test['data_type'] <- 'test'
    
    x_train['target'] <- y_train
    x_test['target'] <- y_test
    
    x <- rbind.data.frame(x_train, x_test)
    x['cancer_type'] <- cancer_type
    x['labeling_type'] <- agg_level
    x['repeat'] <- i
    return(x)
  }
  write.csv(repeats_res, file = paste0(output_path, target_column, ".csv"))
  pb$tick()
}


############################################ BREAKPOINTS RANDOMNESS DATA

output_path <- "data/adhoc/datasets/randomness/"

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


# join target
all_target_data <- q_target %>%
  inner_join(
    hsp_bkpt_target, by=c("chr", "from", "to")
  ) %>%
  select(c("chr", "from", "to", all_target_cols))


# train/test split
train_ratio <- 0.7
n_repeats <- 30

# set progress bar
n_iterations <- length(all_target_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)

for (target_column in all_target_cols){
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  agg_level <- strsplit(target_column, "_")[[1]][2]

  # if it is classification task for "breakpoints" vs "hotspots" then drop rows without any breakpoints
  if (nrow(unique(all_target_data[target_column])) == 3){
    new_all_data <- all_target_data[all_target_data[target_column] != 2, ]
  } else {
    new_all_data <- all_target_data
  }
  
  repeats_res <- foreach(i=seq(1, n_repeats),.combine=rbind) %dopar% {
    
    # train/test split
    splitted_dataset <- get_train_test_split(data=new_all_data, target_col=target_column, 
                                             start_pos_col="from",
                                             chr_col="chr", feature_cols=c("chr", "from", "to"),
                                             train_ratio=train_ratio, 
                                             seed=i)
    x_train <- splitted_dataset[["x_train"]]
    y_train <- splitted_dataset[["y_train"]]
    x_test <- splitted_dataset[["x_test"]] 
    y_test <- splitted_dataset[["y_test"]]  
    
    # save into one dataset
    x_train['data_type'] <- 'train'
    x_test['data_type'] <- 'test'
    
    x_train['target'] <- y_train
    x_test['target'] <- y_test
    
    x <- rbind.data.frame(x_train, x_test)
    x['cancer_type'] <- cancer_type
    x['labeling_type'] <- agg_level
    x['repeat'] <- i
    return(x)
  }
  write.csv(repeats_res, file = paste0(output_path, target_column, ".csv"))
  pb$tick()
}
