script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(reshape2)

library(ggplot2)

library(doParallel)

library(caret)
library(randomForest)
library(pROC)

source("run/features.R")
source("run/tools.R")

n_cores <- 14
registerDoParallel(n_cores)

set.seed(7)


win_len <- 100000

# load data
data <- read.csv(
  paste0("data/datasets/dataset_", format(win_len, scientific = FALSE), ".csv")
)

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

# get new features
all_data_bin <- get_binary_features(data, features_cols)
bin_features <- setdiff(names(all_data_bin), c(features_cols, ss_cols, hsp_cols))

all_data_max <- get_maximum_features(data, features_cols = features_cols)
max_features <- setdiff(names(all_data_max), c(features_cols, ss_cols, hsp_cols))

all_data_high <- get_higher_level_features(data=data, features_cols = features_cols, win_len_upper = 1000000,
                                      path_to_upper_data = "data/datasets/dataset_1000000.csv")
high_features <- setdiff(names(all_data_high), c(features_cols, ss_cols, hsp_cols))

# transform binary features to category
new_feat_names <- c(high_features, max_features, bin_features)
numerical_feats <- c(high_features, grep(x = new_feat_names, pattern = "ratio", value = TRUE))
binary_feats <- setdiff(new_feat_names, numerical_feats)
all_numerical_features <- c(features_cols, numerical_feats)
all_features <- c(features_cols, new_feat_names)

# join all features
all_data <- all_data_bin %>%
  inner_join(
    all_data_max %>% select(c("chr","from", max_features)),
    by=c("chr","from")
    ) %>%
  inner_join(
    all_data_high %>% select(c("chr","from", high_features)),
    by=c("chr","from")
  ) 
rm(all_data_bin)
rm(all_data_max)
rm(all_data_high)
gc()

# binary features to factor
all_data[binary_feats] <- apply(all_data[binary_feats], MARGIN = 2, FUN = function(x) as.factor(as.character(x)))

# drop constant features
feat_info <- data.frame(apply(all_data[all_features], MARGIN = 2, FUN = function(x) {length(unique(x))}))
names(feat_info) <- "counts"
feat_info$feat <- rownames(feat_info)
final_feats <- setdiff(all_features, feat_info[feat_info$counts == 1, 'feat'])
final_feats

all_tissue_spec_feats <- vector()
for (feat in final_feats){
  for (col in tissue_spec_feats){
    all_tissue_spec_feats <- c(all_tissue_spec_feats, grep(x = feat, pattern = col, value = TRUE))  
  }
}
all_conserved_feats <- setdiff(final_feats, all_tissue_spec_feats)

# select same aggregation level to compare with base model
# optimal model parameters will be used

# train_test_split
agg_level <- "_99._"
init_ratio <- 0.5
train_ratio <- 0.7

chrs <-
  c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )
all_data$chr <- as.character(all_data$chr)
all_pos_map <- data.frame()
for (chr in chrs){
  starts <- all_data[all_data$chr == chr, 'from']
  starts <- starts[order(starts)]
  bin_name <- cut(starts, breaks = 4, labels = c("first_bin","second_bin","third_bin","last_bin"),
                  include.lowest = TRUE, right = FALSE)
  pos_map <- data.frame(from = starts, pos_bin = bin_name, chr = chr)
  all_pos_map <- rbind(all_pos_map, pos_map)
}
all_pos_map$chr <- as.character(all_pos_map$chr)
all_data <- all_data %>%
  inner_join(all_pos_map, by=c("chr","from"))

all_data$strata <- paste0(all_data$chr, "_", all_data$pos_bin)

all_target_cols <- grep(x =  hsp_cols, pattern = agg_level, value = TRUE)

all_data$sum_targets <- rowSums(all_data[, all_target_cols])
all_neg_label <- all_data[all_data$sum_targets == 0, ]

set.seed(7)
tr_neg <- createDataPartition(all_neg_label$strata, p = init_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

all_results <- data.frame()

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- all_data[all_data[target_column] == 1, ]
  set.seed(7)
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  
  all_features_cols <- c(all_conserved_feats, grep(x = all_tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  # split on train/test
  n_repeats <- 30
  
  # models parameters
  trCtrl <- trainControl(
    method="none",
    verboseIter = TRUE,
    classProbs=TRUE,
    summaryFunction = twoClassSummary,
    seeds = seq(1, n_repeats))
  # set "seeds" parameters to make results reproducible
  mtryGrid <- expand.grid(mtry = 5)
  
  one_iter_res <- foreach(i=seq(1, n_repeats), .combine = rbind) %dopar% {

    set.seed(i)
    tr_neg <- createDataPartition(neg_data_train$strata, p = train_ratio, list = FALSE)
    data_train <- neg_data_train[tr_neg,]
    data_test  <- neg_data_train[-tr_neg,]
    
    set.seed(i)
    tr_pos <- sample(x = rownames(pos_data_train), size = floor(train_ratio * nrow(pos_data_train)), replace = FALSE)
    te_pos <- setdiff(rownames(pos_data_train), tr_pos)
    data_train <- rbind(data_train, pos_data_train[tr_pos, ])
    data_test <- rbind(data_test, pos_data_train[te_pos, ])
    
    x_train <- data_train[all_features_cols]
    x_test <- data_test[all_features_cols]
    
    y_train <- as.factor(data_train[, target_column])
    levels(y_train) <- c("X0", "X1")
    y_test <- as.factor(data_test[, target_column])
    levels(y_test) <- c("X0", "X1")
    
    # limit outliers
    all_numerical_cols <- intersect(all_features_cols, all_numerical_features)
    data_tr <- limit_outliers(x_train = x_train, x_test = x_test, 
                              features = all_numerical_cols, iqr_multiplicator = 3)
    all_bin_cols <- setdiff(all_features_cols, all_numerical_cols)
    
    x_train <- x_train %>% 
      select(all_bin_cols) %>%
      cbind(data_tr[['train']])
      
    x_test <- x_test %>% 
      select(all_bin_cols) %>%
      cbind(data_tr[['test']])

    
    feat_groups <- list()
    feat_groups[["base"]] <- intersect(features_cols, all_features_cols) 
    feat_groups[["base_with_binary_flags"]] <- intersect(
      c(features_cols, bin_features), 
      all_features_cols)
    feat_groups[["base_with_higher_level_feats"]] <- intersect(
      c(features_cols, high_features), 
      all_features_cols) 
    feat_groups[["base_with_maximums"]] <- intersect(
      c(features_cols, max_features), 
      all_features_cols) 

    feat_iter_res <- data.frame()
    
    # for each feature group
    
    for (j in 1:length(feat_groups)){
      
      features_gr <- names(feat_groups)[j]
      features_nm <- feat_groups[[features_gr]]
      features_nm
      # fit models
      model <- train(
        x = x_train[features_nm], 
        y = y_train,
        preProcess = c("zv", "YeoJohnson", "center", "scale"),
        method = "rf",
        metric="ROC",   
        maximize = TRUE,
        trControl = trCtrl,
        ntree = 1000,
        nodesize = 70,
        maxnodes = 10,
        sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
        tuneGrid = mtryGrid
      )
      
      # ROC AUC
      train_pred <- predict(model, newdata = x_train[features_nm], type = "prob")
      test_pred <- predict(model, newdata = x_test[features_nm], type = "prob")
      train_pred$target <- y_train
      test_pred$target <- y_test

      tr_auc <- auc(response = train_pred$target, predictor = train_pred$X1, levels=c("X0", "X1"), direction = "<") 
      te_auc <- auc(response = test_pred$target, predictor = test_pred$X1, levels=c("X0", "X1"), direction = "<")
      
      out_df <- data.frame(features_group = features_gr, tr_auc = tr_auc, te_auc = te_auc)     
      
      feat_iter_res <- rbind(feat_iter_res, out_df)
    }
    
    feat_iter_res$iter <- i
    
    return(feat_iter_res)
  }
  
  one_iter_res$cancer_type <- cancer_type  
  all_results <- rbind(all_results, one_iter_res)
  
  
}

all_results$diff_auc <- all_results$tr_auc - all_results$te_auc
write.csv(all_results, file = "data/output/feature_engineering.csv")



all_results <- read.csv(file = "data/output/feature_engineering.csv")

all_results_reshaped <- dcast(data = all_results, formula = iter + cancer_type ~ features_group, value.var = "te_auc")

for (col in c( "base_with_binary_flags" , "base_with_higher_level_feats", "base_with_maximums")){
  all_results_reshaped <- all_results_reshaped %>%
    mutate(!!as.name(paste0(col, "_diff")) := !!as.name(col) - base)
}
all_results_reshaped_plot <- melt(data = all_results_reshaped, id.vars = c("cancer_type", "iter"), 
                                  measure.vars = c( "base_with_binary_flags_diff" , "base_with_higher_level_feats_diff", "base_with_maximums_diff"))
ggplot(all_results_reshaped_plot, aes(value, fill=variable)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~cancer_type)+
  geom_vline(xintercept = 0)

mean_diff <- all_results_reshaped %>%
  group_by(cancer_type) %>%
  summarize(
    mean(base_with_binary_flags_diff),
    mean(base_with_higher_level_feats_diff),
    mean(base_with_maximums_diff)
  )









### CHECK BINARY ONLY

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- all_data[all_data[target_column] == 1, ]
  set.seed(7)
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  
  all_features_cols <- c(all_conserved_feats, grep(x = all_tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  # split on train/test
  n_repeats <- 30
  
  # models parameters
  trCtrl <- trainControl(
    method="none",
    verboseIter = TRUE,
    classProbs=TRUE,
    summaryFunction = twoClassSummary,
    seeds = seq(1, n_repeats))
  # set "seeds" parameters to make results reproducible
  mtryGrid <- expand.grid(mtry = 5)
  
  one_iter_res <- foreach(i=seq(1, n_repeats), .combine = rbind) %dopar% {
    
    set.seed(i)
    tr_neg <- createDataPartition(neg_data_train$strata, p = train_ratio, list = FALSE)
    data_train <- neg_data_train[tr_neg,]
    data_test  <- neg_data_train[-tr_neg,]
    
    set.seed(i)
    tr_pos <- sample(x = rownames(pos_data_train), size = floor(train_ratio * nrow(pos_data_train)), replace = FALSE)
    te_pos <- setdiff(rownames(pos_data_train), tr_pos)
    data_train <- rbind(data_train, pos_data_train[tr_pos, ])
    data_test <- rbind(data_test, pos_data_train[te_pos, ])
    
    x_train <- data_train[all_features_cols]
    x_test <- data_test[all_features_cols]
    
    y_train <- as.factor(data_train[, target_column])
    levels(y_train) <- c("X0", "X1")
    y_test <- as.factor(data_test[, target_column])
    levels(y_test) <- c("X0", "X1")
    
    # limit outliers
    all_numerical_cols <- intersect(all_features_cols, all_numerical_features)
    data_tr <- limit_outliers(x_train = x_train, x_test = x_test, 
                              features = all_numerical_cols, iqr_multiplicator = 3)
    all_bin_cols <- setdiff(all_features_cols, all_numerical_cols)
    
    x_train <- x_train %>% 
      select(all_bin_cols) %>%
      cbind(data_tr[['train']])
    
    x_test <- x_test %>% 
      select(all_bin_cols) %>%
      cbind(data_tr[['test']])
    
    
    feat_groups <- list()
    feat_groups[["base"]] <- intersect(features_cols, all_features_cols) 
    feat_groups[["binary_flags"]] <- intersect(bin_features, all_features_cols)
    feat_groups[["binary_with_maximums"]] <- setdiff(c(intersect(bin_features, all_features_cols),
                                               intersect(max_features, all_features_cols)), all_numerical_cols)
    feat_iter_res <- data.frame()
    
    # for each feature group
    
    for (j in 1:length(feat_groups)){
      
      features_gr <- names(feat_groups)[j]
      features_nm <- feat_groups[[features_gr]]

      # fit models
      model <- train(
        x = x_train[features_nm],
        y = y_train,
        preProcess = c("zv"),
        method = "rf",
        metric="ROC",
        maximize = TRUE,
        trControl = trCtrl,
        ntree = 1000,
        nodesize = 70,
        maxnodes = 10,
        sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
        tuneGrid = mtryGrid
      )
      
      # ROC AUC
      train_pred <- predict(model, newdata = x_train[features_nm], type = "prob")
      test_pred <- predict(model, newdata = x_test[features_nm], type = "prob")
      train_pred$target <- y_train
      test_pred$target <- y_test
      
      tr_auc <- auc(response = train_pred$target, predictor = train_pred$X1, levels=c("X0", "X1"), direction = "<") 
      te_auc <- auc(response = test_pred$target, predictor = test_pred$X1, levels=c("X0", "X1"), direction = "<")
      
      out_df <- data.frame(features_group = features_gr, tr_auc = tr_auc, te_auc = te_auc)     
      
      feat_iter_res <- rbind(feat_iter_res, out_df)
    }
    
    feat_iter_res$iter <- i
    
    return(feat_iter_res)
  }
  
  one_iter_res$cancer_type <- cancer_type  
  all_results <- rbind(all_results, one_iter_res)
  
  
}

all_results$diff_auc <- all_results$tr_auc - all_results$te_auc
write.csv(all_results, file = "data/output/feature_engineering_only_binary.csv")


all_results_reshaped <- dcast(data = all_results, formula = iter + cancer_type ~ features_group, value.var = "te_auc")

for (col in c( "binary_flags" , "binary_with_maximums")){
  all_results_reshaped <- all_results_reshaped %>%
    mutate(!!as.name(paste0(col, "_diff")) := !!as.name(col) - base)
}
all_results_reshaped_plot <- melt(data = all_results_reshaped, id.vars = c("cancer_type", "iter"), 
                                  measure.vars = c( "binary_with_maximums_diff", "binary_flags_diff"))
ggplot(all_results_reshaped_plot, aes(value, fill=variable)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~cancer_type)+
  geom_vline(xintercept = 0)

mean_diff <- all_results_reshaped %>%
  group_by(cancer_type) %>%
  summarize(
    mean(binary_with_maximums_diff),
    mean(binary_flags_diff)
  )
