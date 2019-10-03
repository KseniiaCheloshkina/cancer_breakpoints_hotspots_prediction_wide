script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(ggplot2)
library(randomForest)
library(caret)
library(pROC)
library(doParallel)
library(reshape2)
library(ROSE)

n_cores <- 14
registerDoParallel(n_cores)

set.seed(7)


win_len <- 1000000

# load data
data <- read.csv(
  paste0("data/datasets/dataset_", format(win_len, scientific = FALSE), ".csv")
)

hsp_cols <- grep(x = names(data), pattern = "hsp",value = TRUE)

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



# BUILD ONE MODEL
# select agg_level and cancer_type
agg_level <- "_99._"
cancer_type <- 'liver'
target_col <- grep(x = grep(x =  hsp_cols, pattern = agg_level, value = TRUE), pattern = cancer_type, value = TRUE)
all_features_cols <- c(conserved_features, grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))


# create strata for train/test split
chrs <-
  c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )
data$chr <- as.character(data$chr)
all_pos_map <- data.frame()
for (chr in chrs){
  starts <- data[data$chr == chr, 'from']
  starts <- starts[order(starts)]
  bin_name <- cut(starts, breaks = 4, labels = c("first_bin","second_bin","third_bin","last_bin"),
      include.lowest = TRUE, right = FALSE)
  pos_map <- data.frame(from = starts, pos_bin = bin_name, chr = chr)
  all_pos_map <- rbind(all_pos_map, pos_map)
}
all_pos_map$chr <- as.character(all_pos_map$chr)
data <- data %>%
  inner_join(all_pos_map, by=c("chr","from"))

data$strata <- paste0(data$chr, "_", data$pos_bin)

train_ratio <- 0.7

rownames(data) <- NULL
pos_label_data <- data[data[target_col] == 1, ]
neg_label_data <- data[data[target_col] == 0, ]

n_repeats <- 100

res_models <- foreach(i = seq(1, n_repeats), .packages=c('randomForest','pROC', 'caret')) %dopar% { 
  
  set.seed(i)
  
  # train/test data
  tr_neg <- createDataPartition(neg_label_data$strata, p = train_ratio, list = FALSE)
  data_train <- neg_label_data[tr_neg,]
  data_test  <- neg_label_data[-tr_neg,]
  
  tr_pos <- sample(x = rownames(pos_label_data), size = floor(train_ratio * nrow(pos_label_data)), replace = FALSE)
  te_pos <- setdiff(rownames(pos_label_data), tr_pos)
  data_train <- rbind(data_train, pos_label_data[tr_pos, ])
  data_test <- rbind(data_test, pos_label_data[te_pos, ])
  
  x_train <- data_train[all_features_cols]
  x_test <- data_test[all_features_cols]
  
  y_train <- as.factor(data_train[, target_col])
  levels(y_train) <- c("X0", "X1")
  y_test <- as.factor(data_test[,target_col])
  levels(y_test) <- c("X0", "X1")
  
  # fit model
  mtryGrid <- expand.grid(mtry = 6)
  model <- train(
    x = x_train, 
    y = y_train,
    preProcess = c("YeoJohnson", "center", "scale"),
    method = "rf",
    metric="ROC",   
    maximize = TRUE,
    trControl = trainControl(
      method="none",
      number = 5,
      verboseIter = TRUE,
      classProbs=TRUE,
      summaryFunction = twoClassSummary),
    ntree = 50,
    nodesize = 20,
    sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
    tuneGrid = mtryGrid
  )

  train_pred <- predict(model, newdata = x_train, type = "prob")
  test_pred <- predict(model, newdata = x_test, type = "prob")
  
  tr_auc <- auc(y_train, train_pred[, 2]) 
  te_auc <- auc(y_test, test_pred[, 2]) 

  return(list(as.numeric(tr_auc), as.numeric(te_auc)))
} 

all_tr_auc <- unlist(lapply(res_models, function(x) x[[1]]))
all_te_auc <- unlist(lapply(res_models, function(x) x[[2]]))

roc_auc_data <- data.frame(train_auc = all_tr_auc, test_auc = all_te_auc)
roc_auc_data_m <- melt(roc_auc_data)

ggplot(roc_auc_data_m, aes(value, fill=variable)) + 
  geom_density(alpha = 0.5)

ggplot(roc_auc_data, aes(x=train_auc, y=test_auc)) + 
  geom_point()





####################################### MODEL SELECTION FOR RANDOM FOREST

agg_level <- "_99._"
init_ratio <- 0.5
train_ratio <- 0.7

chrs <-
  c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )
data$chr <- as.character(data$chr)
all_pos_map <- data.frame()
for (chr in chrs){
  starts <- data[data$chr == chr, 'from']
  starts <- starts[order(starts)]
  bin_name <- cut(starts, breaks = 4, labels = c("first_bin","second_bin","third_bin","last_bin"),
                  include.lowest = TRUE, right = FALSE)
  pos_map <- data.frame(from = starts, pos_bin = bin_name, chr = chr)
  all_pos_map <- rbind(all_pos_map, pos_map)
}
all_pos_map$chr <- as.character(all_pos_map$chr)
data <- data %>%
  inner_join(all_pos_map, by=c("chr","from"))

data$strata <- paste0(data$chr, "_", data$pos_bin)

all_target_cols <- grep(x =  hsp_cols, pattern = agg_level, value = TRUE)

data$sum_targets <- rowSums(data[, all_target_cols])
all_neg_label <- data[data$sum_targets == 0, ]


tr_neg <- createDataPartition(all_neg_label$strata, p = init_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

all_results <- data.frame()

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- data[data[target_column] == 1, ]
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]

  cancer_type <- strsplit(target_column, "_")[[1]][3]
  all_features_cols <- c(conserved_features, grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  # split on train/test
  
  n_repeats <- 10
  
  par_df <- expand.grid(mtry = c(3, 4, 5, 6), 
                        ntree = c(500, 1000), 
                        nodesize = c(30, 50, 70),
                        maxnodes = c(3, 5, 8, 10, 12)
                        )
  
  for (i in seq(1, n_repeats)){
    
    set.seed(i)
    
    tr_neg <- createDataPartition(neg_data_train$strata, p = train_ratio, list = FALSE)
    data_train <- neg_data_train[tr_neg,]
    data_test  <- neg_data_train[-tr_neg,]
    
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
    
    # fit models
    
    res_iter <- foreach(i=seq(1, nrow(par_df)), .combine = rbind) %dopar% {
      
      mtry <- par_df[i, ]$mtry
      ntree <- par_df[i, ]$ntree
      nodesize <- par_df[i, ]$nodesize
      maxnodes <- par_df[i, ]$maxnodes
      
      mtryGrid <- expand.grid(mtry = mtry)
      model <- train(
        x = x_train, 
        y = y_train,
        preProcess = c("YeoJohnson", "center", "scale"),
        method = "rf",
        metric="ROC",   
        maximize = TRUE,
        trControl = trainControl(
          method="none",
          verboseIter = TRUE,
          classProbs=TRUE,
          summaryFunction = twoClassSummary),
        ntree = ntree,
        nodesize = nodesize,
        maxnodes = maxnodes,
        sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
        tuneGrid = mtryGrid
      )
      
      train_pred <- predict(model, newdata = x_train, type = "prob")
      test_pred <- predict(model, newdata = x_test, type = "prob")
      
      tr_auc <- auc(y_train, train_pred[, 2]) 
      te_auc <- auc(y_test, test_pred[, 2]) 
      
      out_df <- data.frame(mtry = mtry, ntree = ntree, nodesize = nodesize, maxnodes = maxnodes, tr_auc = tr_auc, te_auc = te_auc, iter = i)
      return(out_df)
    }
    
    res_iter$iter <- i
    res_iter$cancer_type <- cancer_type
    all_results <- rbind(all_results, res_iter)
    
  }

}

all_results$diff_auc <- all_results$tr_auc - all_results$te_auc
write.csv(all_results, file = "data/hp_rf.csv")

all_results <- read.csv("data/hp_rf.csv")
res_all <- all_results %>%
  group_by(cancer_type) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q20 = quantile(te_auc, 0.2),
    q30 = quantile(te_auc, 0.3),
    sd = sd(te_auc)
  ) %>%
  arrange(cancer_type)

pred_cancers_data <- all_results %>%
  filter(cancer_type %in% c("breast", "ovary", "prostate", "brain")) 

ggplot(pred_cancers_data, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)

ggplot(all_results, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)

all_results %>%
  filter(cancer_type %in% c('brain', 'breast', 'ovary','prostate')) %>%
  group_by(ntree) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q30 = quantile(te_auc, 0.3),
    q70 = quantile(te_auc, 0.7),
    max_te = max(te_auc)
  ) %>%
  arrange(min_te_auc)

gb_pars <- pred_cancers_data %>% 
  group_by(ntree, nodesize, maxnodes, mtry, cancer_type) %>%
  summarize(
    mean_d = mean(diff_auc),
    med = median(diff_auc),
    sd = sd(te_auc)
  ) %>%
  group_by(cancer_type) %>%
  mutate(
    mean_rank = dense_rank(mean_d),
    med_rank = dense_rank(med),
    sd_rank = ntile(sd, 5)
  ) 
gb_pars <- gb_pars %>%
  group_by(ntree, nodesize, maxnodes, mtry) %>%
  summarize(
    mean_rank = mean(mean_rank),
    med_rank = mean(med_rank),
    sd_rank = mean(sd_rank)
  ) %>%
  mutate(summed_rank = mean_rank + med_rank + sd_rank)

gb_pars %>%
  ungroup() %>%
  filter(summed_rank == min(summed_rank))

opt_pars <- all_results %>%
  filter(
    nodesize == 50,
    ntree == 500,
    maxnodes == 10,
    mtry == 6
  )

res_opt <- opt_pars %>%
  group_by(cancer_type) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q20 = quantile(te_auc, 0.2),
    med_te_auc = median(te_auc),
    mean_te_auc = mean(te_auc),
    q80 = quantile(te_auc, 0.8),
    q90 = quantile(te_auc, 0.9),
    max_te_auc = max(te_auc),
    sd = sd(te_auc)
  ) %>%
  arrange(cancer_type)

ggplot(opt_pars, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)



######## REMOVE CORRELATED FEATURES WITH PLS AND AGAIN TUNE RANDOM FOREST
library(pls)
# https://www.kaggle.com/sasali/notebook83cb3b1069#PLS---RF(customized-caret-train-function)

PlsRf <- list(label = "PLS-RF", type = "Classification",
              library = c("pls", "randomForest"),
              loop = NULL,
              parameters = data.frame(parameter = c('ncomp', 'ntree', 'mtry', 'nodesize','maxnodes'),
                                      class = c("numeric", 'numeric', 'numeric', 'numeric', 'numeric'),
                                      label = c('#Components', '#Trees',
                                                '#Randomly Selected Predictors',  
                                                '#Maximum size of terminal nodes',
                                                '#Maximum number of terminal nodes')),
              grid = function(x, y, len = NULL) {
                grid <- expand.grid(ncomp = seq(1, min(ncol(x) - 1, len), by = 1),
                                    mtry = 1:len)
                grid <- subset(grid, mtry <= ncomp)
              },
              
              fit = function(x, y, wts, param, lev, last, classProbs, ...) { 
                # The probability model is fit based on the value of classProbs. 
                # This value is determined by the value given in trainControl
                pre <- plsda(x, y, ncomp = param$ncomp)
                scores <- pls:::predict.mvr(pre, x, type = "scores")
                mod <- randomForest(scores, y, mtry = param$mtry, ntree=param$ntree, 
                                    nodesize=param$nodesize, maxnodes=param$maxnodes)
                mod$projection <- pre$projection
                mod },
              
              predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {       
                scores <- as.matrix(newdata) %*% modelFit$projection
                predict(modelFit, scores)
              },
              
              prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {       
                scores <- as.matrix(newdata) %*% modelFit$projection
                predict(modelFit, scores, type="prob")
              },
              
              varImp = NULL,
              predictors = function(x, ...) rownames(x$projection),
              levels = function(x) x$obsLevels,
              sort = function(x) x[order(x[,1]),]   
)

agg_level <- "_99._"
init_ratio <- 0.5
train_ratio <- 0.7

chrs <-
  c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )
data$chr <- as.character(data$chr)
all_pos_map <- data.frame()
for (chr in chrs){
  starts <- data[data$chr == chr, 'from']
  starts <- starts[order(starts)]
  bin_name <- cut(starts, breaks = 4, labels = c("first_bin","second_bin","third_bin","last_bin"),
                  include.lowest = TRUE, right = FALSE)
  pos_map <- data.frame(from = starts, pos_bin = bin_name, chr = chr)
  all_pos_map <- rbind(all_pos_map, pos_map)
}
all_pos_map$chr <- as.character(all_pos_map$chr)
data <- data %>%
  inner_join(all_pos_map, by=c("chr","from"))

data$strata <- paste0(data$chr, "_", data$pos_bin)

all_target_cols <- grep(x =  hsp_cols, pattern = agg_level, value = TRUE)

data$sum_targets <- rowSums(data[, all_target_cols])
all_neg_label <- data[data$sum_targets == 0, ]


tr_neg <- createDataPartition(all_neg_label$strata, p = init_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

all_results <- data.frame()

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- data[data[target_column] == 1, ]
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  all_features_cols <- c(conserved_features, grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  # split on train/test
  
  n_repeats <- 10
  
  par_df <- expand.grid(mtry = c(2, 3, 4), 
                        ntree = c(500, 1000), 
                        nodesize = c(50, 70),
                        maxnodes = c(4, 6, 8),
                        ncomp = c(5, 10, 15)
  )
  
  for (i in seq(1, n_repeats)){
    
    set.seed(i)
    
    tr_neg <- createDataPartition(neg_data_train$strata, p = train_ratio, list = FALSE)
    data_train <- neg_data_train[tr_neg,]
    data_test  <- neg_data_train[-tr_neg,]
    
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
    
    # fit models
    
    res_iter <- foreach(i=seq(1, nrow(par_df)), .combine = rbind) %dopar% {
      
      mtry <- par_df[i, ]$mtry
      ntree <- par_df[i, ]$ntree
      nodesize <- par_df[i, ]$nodesize
      maxnodes <- par_df[i, ]$maxnodes
      ncomp <- par_df[i, ]$ncomp
      
      mtryGrid <- expand.grid(
        mtry = mtry,
        ntree = ntree,
        nodesize = nodesize,
        maxnodes = maxnodes,
        ncomp = ncomp
        )
      
      model <- train(
        x = x_train, 
        y = y_train,
        preProcess = c("YeoJohnson", "center", "scale"),
        method = PlsRf,
        metric="ROC",   
        maximize = TRUE,
        trControl = trainControl(
          method="none",
          verboseIter = TRUE,
          classProbs=TRUE,
          summaryFunction = twoClassSummary),
        sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
        tuneGrid = mtryGrid
      )
      
      train_pred <- predict(model, newdata = x_train, type = "prob")
      test_pred <- predict(model, newdata = x_test, type = "prob")
      
      tr_auc <- auc(y_train, train_pred[, 2]) 
      te_auc <- auc(y_test, test_pred[, 2]) 
      
      out_df <- data.frame(mtry = mtry, ntree = ntree, nodesize = nodesize, maxnodes = maxnodes, ncomp = ncomp, tr_auc = tr_auc, te_auc = te_auc, iter = i)
      return(out_df)
    }
    
    res_iter$iter <- i
    res_iter$cancer_type <- cancer_type
    all_results <- rbind(all_results, res_iter)
    
  }
  
}

all_results$diff_auc <- all_results$tr_auc - all_results$te_auc
write.csv(all_results, file = "data/hp_rf_pls.csv")


all_results <- read.csv("data/hp_rf_pls.csv")
res_all <- all_results %>%
  group_by(cancer_type) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q20 = quantile(te_auc, 0.2),
    q30 = quantile(te_auc, 0.3),
    sd = sd(te_auc)
  ) %>%
  arrange(cancer_type)

pred_cancers_data <- all_results %>%
  filter(cancer_type %in% c("breast", "ovary", "prostate", "brain")) 

ggplot(pred_cancers_data, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)

ggplot(all_results, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)

all_results %>%
  filter(cancer_type %in% c('brain', 'breast', 'ovary','prostate')) %>%
  group_by(ntree) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q30 = quantile(te_auc, 0.3),
    q70 = quantile(te_auc, 0.7),
    max_te = max(te_auc)
  ) %>%
  arrange(min_te_auc)

gb_pars <- pred_cancers_data %>% 
  group_by(ntree, nodesize, maxnodes, mtry, ncomp, cancer_type) %>%
  summarize(
    mean_d = mean(diff_auc),
    med = median(diff_auc),
    sd = sd(te_auc)
  ) %>%
  group_by(cancer_type) %>%
  mutate(
    mean_rank = dense_rank(mean_d),
    med_rank = dense_rank(med),
    sd_rank = ntile(sd, 5)
  ) 
gb_pars <- gb_pars %>%
  group_by(ntree, nodesize, maxnodes, mtry, ncomp) %>%
  summarize(
    mean_rank = mean(mean_rank),
    med_rank = mean(med_rank),
    sd_rank = mean(sd_rank)
  ) %>%
  mutate(summed_rank = mean_rank + med_rank + sd_rank)

gb_pars %>%
  ungroup() %>%
  filter(summed_rank == min(summed_rank))

opt_pars <- all_results %>%
  filter(
    nodesize == 50,
    ntree == 1000,
    maxnodes == 8,
    mtry == 3,
    ncomp == 15
  )

res_opt <- opt_pars %>%
  group_by(cancer_type) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q20 = quantile(te_auc, 0.2),
    med_te_auc = median(te_auc),
    mean_te_auc = mean(te_auc),
    q80 = quantile(te_auc, 0.8),
    q90 = quantile(te_auc, 0.9),
    max_te_auc = max(te_auc),
    sd = sd(te_auc)
  ) %>%
  arrange(cancer_type)

ggplot(opt_pars, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)



######## REMOVE CORRELATED FEATURES WITH PCA AND AGAIN TUNE RANDOM FOREST

agg_level <- "_99._"
init_ratio <- 0.5
train_ratio <- 0.7

chrs <-
  c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )
data$chr <- as.character(data$chr)
all_pos_map <- data.frame()
for (chr in chrs){
  starts <- data[data$chr == chr, 'from']
  starts <- starts[order(starts)]
  bin_name <- cut(starts, breaks = 4, labels = c("first_bin","second_bin","third_bin","last_bin"),
                  include.lowest = TRUE, right = FALSE)
  pos_map <- data.frame(from = starts, pos_bin = bin_name, chr = chr)
  all_pos_map <- rbind(all_pos_map, pos_map)
}
all_pos_map$chr <- as.character(all_pos_map$chr)
data <- data %>%
  inner_join(all_pos_map, by=c("chr","from"))

data$strata <- paste0(data$chr, "_", data$pos_bin)

all_target_cols <- grep(x =  hsp_cols, pattern = agg_level, value = TRUE)

data$sum_targets <- rowSums(data[, all_target_cols])
all_neg_label <- data[data$sum_targets == 0, ]


tr_neg <- createDataPartition(all_neg_label$strata, p = init_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

all_results <- data.frame()

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- data[data[target_column] == 1, ]
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  all_features_cols <- c(conserved_features, grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  # split on train/test
  
  n_repeats <- 10
  
  par_df <- expand.grid(mtry = c(2, 3, 4), 
                        ntree = c(500, 1000), 
                        nodesize = c(50, 70),
                        maxnodes = c(8, 10, 12),
                        pcaComp = c(5, 10, 15)
  )
  
  for (i in seq(1, n_repeats)){
    
    set.seed(i)
    
    tr_neg <- createDataPartition(neg_data_train$strata, p = train_ratio, list = FALSE)
    data_train <- neg_data_train[tr_neg,]
    data_test  <- neg_data_train[-tr_neg,]
    
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
    
    # fit models
    
    res_iter <- foreach(i=seq(1, nrow(par_df)), .combine = rbind) %dopar% {
      
      mtry <- par_df[i, ]$mtry
      ntree <- par_df[i, ]$ntree
      nodesize <- par_df[i, ]$nodesize
      maxnodes <- par_df[i, ]$maxnodes
      pcaComp <- par_df[i, ]$pcaComp
      
      mtryGrid <- expand.grid(mtry = mtry)
      model <- train(
        x = x_train, 
        y = y_train,
        preProcess = c("YeoJohnson", "center", "scale", "pca"),
        method = "rf",
        metric="ROC",   
        maximize = TRUE,
        trControl = trainControl(
          method="none",
          verboseIter = TRUE,
          classProbs=TRUE,
          summaryFunction = twoClassSummary,
          preProcOptions = list(pcaComp = pcaComp)
          ),
        ntree = ntree,
        nodesize = nodesize,
        maxnodes = maxnodes,
        sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
        tuneGrid = mtryGrid
      )
      
      train_pred <- predict(model, newdata = x_train, type = "prob")
      test_pred <- predict(model, newdata = x_test, type = "prob")
      
      tr_auc <- auc(y_train, train_pred[, 2]) 
      te_auc <- auc(y_test, test_pred[, 2]) 
      
      out_df <- data.frame(mtry = mtry, ntree = ntree, nodesize = nodesize, maxnodes = maxnodes, pcaComp = pcaComp,
                           tr_auc = tr_auc, te_auc = te_auc, iter = i)
      return(out_df)
    }
    
    res_iter$iter <- i
    res_iter$cancer_type <- cancer_type
    all_results <- rbind(all_results, res_iter)
    
  }
  
}
all_results$diff_auc <- all_results$tr_auc - all_results$te_auc
write.csv(all_results, file = "data/hp_rf_pca.csv")

all_results <- read.csv("data/hp_rf_pca.csv")
res_all <- all_results %>%
  group_by(cancer_type) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q20 = quantile(te_auc, 0.2),
    q30 = quantile(te_auc, 0.3),
    sd = sd(te_auc)
  ) %>%
  arrange(cancer_type)

pred_cancers_data <- all_results %>%
  filter(cancer_type %in% c("breast", "ovary", "prostate", "brain")) 

ggplot(pred_cancers_data, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)

ggplot(all_results, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)

all_results %>%
  filter(cancer_type %in% c('brain', 'breast', 'ovary','prostate')) %>%
  group_by(ntree) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q30 = quantile(te_auc, 0.3),
    q70 = quantile(te_auc, 0.7),
    max_te = max(te_auc)
  ) %>%
  arrange(min_te_auc)

gb_pars <- pred_cancers_data %>% 
  group_by(ntree, nodesize, maxnodes, mtry, pcaComp, cancer_type) %>%
  summarize(
    mean_d = mean(diff_auc),
    med = median(diff_auc),
    sd = sd(te_auc)
  ) %>%
  group_by(cancer_type) %>%
  mutate(
    mean_rank = dense_rank(mean_d),
    med_rank = dense_rank(med),
    sd_rank = ntile(sd, 5)
  ) 
gb_pars <- gb_pars %>%
  group_by(ntree, nodesize, maxnodes, mtry, pcaComp) %>%
  summarize(
    mean_rank = mean(mean_rank),
    med_rank = mean(med_rank),
    sd_rank = mean(sd_rank)
  ) %>%
  mutate(summed_rank = mean_rank + med_rank + sd_rank)

gb_pars %>%
  ungroup() %>%
  filter(summed_rank == min(summed_rank))

opt_pars <- all_results %>%
  filter(
    nodesize == 50,
    ntree == 500,
    maxnodes == 8,
    mtry == 3,
    pcaComp == 5
  )

res_opt <- opt_pars %>%
  group_by(cancer_type) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q20 = quantile(te_auc, 0.2),
    med_te_auc = median(te_auc),
    mean_te_auc = mean(te_auc),
    q80 = quantile(te_auc, 0.8),
    q90 = quantile(te_auc, 0.9),
    max_te_auc = max(te_auc),
    sd = sd(te_auc)
  ) %>%
  arrange(cancer_type)

ggplot(opt_pars, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)



########################### ADABOOST
library(ada)
# discrete adaboost with classification trees, exponential loss 

agg_level <- "_99._"
init_ratio <- 0.5
train_ratio <- 0.7

chrs <-
  c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )
data$chr <- as.character(data$chr)
all_pos_map <- data.frame()
for (chr in chrs){
  starts <- data[data$chr == chr, 'from']
  starts <- starts[order(starts)]
  bin_name <- cut(starts, breaks = 4, labels = c("first_bin","second_bin","third_bin","last_bin"),
                  include.lowest = TRUE, right = FALSE)
  pos_map <- data.frame(from = starts, pos_bin = bin_name, chr = chr)
  all_pos_map <- rbind(all_pos_map, pos_map)
}
all_pos_map$chr <- as.character(all_pos_map$chr)
data <- data %>%
  inner_join(all_pos_map, by=c("chr","from"))

data$strata <- paste0(data$chr, "_", data$pos_bin)

all_target_cols <- grep(x =  hsp_cols, pattern = agg_level, value = TRUE)

data$sum_targets <- rowSums(data[, all_target_cols])
all_neg_label <- data[data$sum_targets == 0, ]


tr_neg <- createDataPartition(all_neg_label$strata, p = init_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

all_results <- data.frame()

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- data[data[target_column] == 1, ]
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  all_features_cols <- c(conserved_features, grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  # split on train/test
  
  n_repeats <- 10
  
  par_df <- expand.grid(nu = c(1, 0.1, 2), 
                        iter = c(10, 50, 100, 500), 
                        maxdepth = c(2, 3, 4)
  )
  
  for (i in seq(1, n_repeats)){
    
    set.seed(i)
    
    tr_neg <- createDataPartition(neg_data_train$strata, p = train_ratio, list = FALSE)
    data_train <- neg_data_train[tr_neg,]
    data_test  <- neg_data_train[-tr_neg,]
    
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
    
    # fit models
    
    res_iter <- foreach(i=seq(1, nrow(par_df)), .combine = rbind) %dopar% {

      iter <- par_df[i, ]$iter
      maxdepth <- par_df[i, ]$maxdepth
      nu <- par_df[i, ]$nu
      
      iterGrid <- expand.grid(iter = iter, maxdepth = maxdepth, nu = nu)
      model <- train(
        x = x_train, 
        y = y_train,
        preProcess = c("YeoJohnson", "center", "scale"),
        method = "ada",
        metric="ROC",   
        maximize = TRUE,
        trControl = trainControl(
          method="none",
          verboseIter = TRUE,
          classProbs=TRUE,
          summaryFunction = twoClassSummary),
        tuneGrid = iterGrid
      )

      train_pred <- predict(model, newdata = x_train, type = "prob")
      test_pred <- predict(model, newdata = x_test, type = "prob")
      
      tr_auc <- auc(y_train, train_pred[, 2]) 
      te_auc <- auc(y_test, test_pred[, 2]) 
      
      out_df <- data.frame(ntrees = iter, maxdepth = maxdepth, nu = nu, tr_auc = tr_auc, te_auc = te_auc, iter = i)
      return(out_df)
    }
    
    res_iter$iter <- i
    res_iter$cancer_type <- cancer_type
    all_results <- rbind(all_results, res_iter)
    
  }
  
}

all_results$diff_auc <- all_results$tr_auc - all_results$te_auc
write.csv(all_results, file = "data/hp_ada.csv")

all_results <- read.csv("data/hp_ada.csv")
res_all <- all_results %>%
  group_by(cancer_type) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q20 = quantile(te_auc, 0.2),
    q30 = quantile(te_auc, 0.3),
    sd = sd(te_auc)
  ) %>%
  arrange(cancer_type)

pred_cancers_data <- all_results %>%
  filter(cancer_type %in% c("breast", "ovary", "prostate", "brain")) 

ggplot(pred_cancers_data, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)

ggplot(all_results, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)

all_results %>%
  filter(cancer_type %in% c('brain', 'breast', 'ovary','prostate')) %>%
  group_by(ntrees) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q30 = quantile(te_auc, 0.3),
    q70 = quantile(te_auc, 0.7),
    max_te = max(te_auc)
  ) %>%
  arrange(min_te_auc)

gb_pars <- pred_cancers_data %>% 
  group_by(ntrees, nu, maxdepth, cancer_type) %>%
  summarize(
    mean_d = mean(diff_auc),
    med = median(diff_auc),
    sd = sd(te_auc)
  ) %>%
  group_by(cancer_type) %>%
  mutate(
    mean_rank = dense_rank(mean_d),
    med_rank = dense_rank(med),
    sd_rank = ntile(sd, 5)
  ) 

gb_pars <- gb_pars %>%
  group_by(ntrees, nu, maxdepth) %>%
  summarize(
    mean_rank = mean(mean_rank),
    med_rank = mean(med_rank),
    sd_rank = mean(sd_rank)
  ) %>%
  mutate(summed_rank = mean_rank + med_rank + sd_rank)

gb_pars %>%
  ungroup() %>%
  filter(summed_rank == min(summed_rank))

opt_pars <- all_results %>%
  filter(
    nu == 0.1,
    ntrees == 50,
    maxdepth == 2
  )

res_opt <- opt_pars %>%
  group_by(cancer_type) %>%
  summarize(
    min_te_auc = min(te_auc),
    q10 = quantile(te_auc, 0.1),
    q20 = quantile(te_auc, 0.2),
    med_te_auc = median(te_auc),
    mean_te_auc = mean(te_auc),
    q80 = quantile(te_auc, 0.8),
    q90 = quantile(te_auc, 0.9),
    max_te_auc = max(te_auc),
    sd = sd(te_auc)
  ) %>%
  arrange(cancer_type)

ggplot(opt_pars, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)



############################### CHECK NEW FEATURES
source("run/features.R")

data <- read.csv("data/datasets/dataset_10000.csv")

hsp_cols <- grep(x = names(data), pattern = "hsp", value = TRUE)
ss_cols <- c ("chr", "from", "to")
features_cols <- setdiff(names(data), c(ss_cols, hsp_cols))

all_data <- get_binary_features(data[features_cols])
# binary as factors???????????????????
all_data <- get_higher_level_features(data=data, features_cols = features_cols, win_len_upper = 100000, 
                                      path_to_upper_data = "data/datasets/dataset_100000.csv")
















############ APPLY BEST SELECTED MODELS ITERATIVELY ON ALL DATA AS USUAL - GET CLASSIFICATION BENCHMARK






############ PU LEARNING
# https://phys.org/news/2018-11-smarter-aimachine-negative.html
# https://www.researchgate.net/publication/4348735_Learning_from_Positive_and_Unlabeled_Examples_A_Survey
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2248178/
# https://link.springer.com/chapter/10.1007/978-3-319-14717-8_45
# http://members.cbio.mines-paristech.fr/~jvert/svn/bibli/local/Pelckmans2009Transductively.pdf
# https://www.comp.nus.edu.sg/~leews/publications/noisyicml.pdf

# https://medium.com/neuralspace/bayesian-neural-network-series-post-1-need-for-bayesian-networks-e209e66b70b2

# https://stackoverflow.com/questions/20150525/stratified-sampling-doesnt-seem-to-change-randomforest-results/20151341#20151341
# https://www.r-bloggers.com/handling-class-imbalance-with-r-and-caret-an-introduction/

# https://topepo.github.io/caret/adaptive-resampling.html
# https://topepo.github.io/caret/feature-selection-using-simulated-annealing.html
# https://topepo.github.io/caret/train-models-by-tag.html#bayesian-model
# https://topepo.github.io/caret/train-models-by-tag.html#oblique-tree
# https://topepo.github.io/caret/models-clustered-by-tag-similarity.html



# ROSE (BAD TRY - GENERATED EXAMPLES ARE TOO DIFFERENT FROM NEGATIVIES AND SIMILAR TO ORIGINAL)
agg_level <- "_99._"
neg_train_ratio <- 0.5
pos_train_ratio <- 0.7

shuffle_train_ratio <- 0.7

chrs <-
  c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )
data$chr <- as.character(data$chr)
all_pos_map <- data.frame()
for (chr in chrs){
  starts <- data[data$chr == chr, 'from']
  starts <- starts[order(starts)]
  bin_name <- cut(starts, breaks = 4, labels = c("first_bin","second_bin","third_bin","last_bin"),
                  include.lowest = TRUE, right = FALSE)
  pos_map <- data.frame(from = starts, pos_bin = bin_name, chr = chr)
  all_pos_map <- rbind(all_pos_map, pos_map)
}
all_pos_map$chr <- as.character(all_pos_map$chr)
data <- data %>%
  inner_join(all_pos_map, by=c("chr","from"))

data$strata <- paste0(data$chr, "_", data$pos_bin)

all_target_cols <- grep(x =  hsp_cols, pattern = agg_level, value = TRUE)

data$sum_targets <- rowSums(data[, all_target_cols])
all_neg_label <- data[data$sum_targets == 0, ]


tr_neg <- createDataPartition(all_neg_label$strata, p = neg_train_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

target_column <- all_target_cols[1]

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- data[data[target_column] == 1, ]
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(pos_train_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]
  all_data_train <- rbind(neg_data_train, pos_data_train)
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  all_features_cols <- c(conserved_features, grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  rose_train <- ROSE(
    formula = as.formula(paste0(target_column, " ~ ", paste0(all_features_cols, collapse = " + "))), 
    data = all_data_train,
    N = nrow(pos_data_train),
    p = 0.9,
    seed = 95
  )$data 
  
  cols <- c(target_column, all_features_cols,"strata")
  rose_train$strata <- "0"
  synth_data <- rbind(rose_train[, cols], neg_data_train[, cols])
  synth_data$strata
  
  # split synthetized data on train/test
  rownames(synth_data) <- NULL
  
  pos_label_data <- synth_data[synth_data[target_column] == 1, ]
  neg_label_data <- synth_data[(synth_data[target_column] == 0) & (synth_data$strata != "0"),  ]
  
  n_repeats <- 10
  #########################
  tr_neg <- createDataPartition(neg_label_data$strata, p = shuffle_train_ratio, list = FALSE)
  data_train <- neg_label_data[tr_neg,]
  data_test  <- neg_label_data[-tr_neg,]
  
  tr_pos <- sample(x = rownames(pos_label_data), size = floor(shuffle_train_ratio * nrow(pos_label_data)), replace = FALSE)
  te_pos <- setdiff(rownames(pos_label_data), tr_pos)
  data_train <- rbind(data_train, pos_label_data[tr_pos, ])
  data_test <- rbind(data_test, pos_label_data[te_pos, ])
  
  x_train <- data_train[all_features_cols]
  x_test <- data_test[all_features_cols]
  
  y_train <- as.factor(data_train[, target_column])
  levels(y_train) <- c("X0", "X1")
  y_test <- as.factor(data_test[, target_column])
  levels(y_test) <- c("X0", "X1")
  
  # fit models
  par_df <- expand.grid(mtry = c(3, 6), 
                        ntree = c(10, 50, 100),
                        nodesize = c(10, 20, 30))
  
  res_iter <- foreach(i=seq(1, nrow(par_df)), .combine = rbind) %dopar% {
    
    mtry <- par_df[i, ]$mtry
    ntree <- par_df[i, ]$ntree
    nodesize <- par_df[i, ]$nodesize
    
    mtryGrid <- expand.grid(mtry = mtry)
    model <- train(
      x = x_train, 
      y = y_train,
      preProcess = c("YeoJohnson", "center", "scale"),
      method = "rf",
      metric="ROC",   
      maximize = TRUE,
      trControl = trainControl(
        method="none",
        verboseIter = TRUE,
        classProbs=TRUE,
        summaryFunction = twoClassSummary),
      ntree = ntree,
      nodesize = nodesize,
      sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
      tuneGrid = mtryGrid
    )
    
    train_pred <- predict(model, newdata = x_train, type = "prob")
    test_pred <- predict(model, newdata = x_test, type = "prob")
    
    tr_auc <- auc(y_train, train_pred[, 2]) 
    te_auc <- auc(y_test, test_pred[, 2]) 
    
    out_df <- data.frame(mtry = mtry, ntree = ntree, nodesize = nodesize, tr_auc = tr_auc, te_auc = te_auc, iter = i)
    return(out_df)
  }
}
