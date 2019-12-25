script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(reshape2)

library(ggplot2)
library(progress)
library(doParallel)

library(caret)
library(randomForest)
library(pls)
library(ada)
library(pROC)

source("run/tools.R")

n_cores <- 15
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



################################## MODEL SELECTION FOR RANDOM FOREST WITHOUT LIMITING OUTLIERS
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
  bin_name <- cut(starts, breaks = 4, labels = c("first_bin", "second_bin", "third_bin", 
                                                 "last_bin"),
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

set.seed(7)
tr_neg <- createDataPartition(all_neg_label$strata, p = init_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

all_results <- data.frame()

# set progress bar
n_iterations <- length(all_target_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- data[data[target_column] == 1, ]
  set.seed(7)
  tr_pos <- sample(rownames(cancer_pos_label_data), 
                   size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]

  cancer_type <- strsplit(target_column, "_")[[1]][3]
  all_features_cols <- c(conserved_features, 
                         grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  par_df <- expand.grid(mtry = c(3, 4, 5, 6), 
                        ntree = c(500, 1000), 
                        nodesize = c(30, 50, 70),
                        maxnodes = c(3, 5, 8, 10, 12)
                        )
  # split on train/test
  n_repeats <- 30
  
  for (i in seq(1, n_repeats)){
    print(i)
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
    # data_tr <- limit_outliers(x_train = x_train, x_test = x_test, 
    #                           features = all_features_cols, iqr_multiplicator = 3)
    # x_train <- data_tr[['train']]
    # x_test <- data_tr[['test']]
    
    # fit models
    trCtrl <- trainControl(
      method="none",
      verboseIter = TRUE,
      classProbs=TRUE,
      summaryFunction = twoClassSummary,
      seeds = seq(1, nrow(par_df)))
    # set "seeds" parameters to make results reproducible
    
    res_iter <- foreach(j=seq(1, nrow(par_df)), .combine = rbind) %dopar% {
      
      mtry <- par_df[j, ]$mtry
      ntree <- par_df[j, ]$ntree
      nodesize <- par_df[j, ]$nodesize
      maxnodes <- par_df[j, ]$maxnodes
      
      mtryGrid <- expand.grid(mtry = mtry)
      model <- train(
        x = x_train, 
        y = y_train,
        preProcess = c("zv", "YeoJohnson", "center", "scale"),
        method = "rf",
        metric="ROC",   
        maximize = TRUE,
        trControl = trCtrl,
        ntree = ntree,
        nodesize = nodesize,
        maxnodes = maxnodes,
        sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
        tuneGrid = mtryGrid
      )

      train_pred <- predict(model, newdata = x_train, type = "prob")
      test_pred <- predict(model, newdata = x_test, type = "prob")
      
      tr_auc <- auc(response = y_train, predictor = train_pred[, 2], levels=c("X0", "X1"),  direction = "<") 
      te_auc <- auc(response = y_test, predictor = test_pred[, 2], levels=c("X0", "X1"),  direction = "<") 

      out_df <- data.frame(
        mtry = mtry, ntree = ntree, nodesize = nodesize, maxnodes = maxnodes,
        tr_auc = tr_auc, te_auc = te_auc, iter = j)
      
      return(out_df)
    }
    
    res_iter$iter <- i
    res_iter$cancer_type <- cancer_type
    all_results <- rbind(all_results, res_iter)
    write.csv(all_results, file = "data/output/model_selection/hp_rf.csv")
    
  }
  pb$tick()
}






####################################### MODEL SELECTION FOR RANDOM FOREST WITH LIMITING OUTLIERS

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

set.seed(7)
tr_neg <- createDataPartition(all_neg_label$strata, p = init_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

all_results <- data.frame()


# all_target_cols <- all_target_cols[all_target_cols != "hsp_99._blood"]

# set progress bar
n_iterations <- length(all_target_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- data[data[target_column] == 1, ]
  set.seed(7)
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  all_features_cols <- c(conserved_features, grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  par_df <- expand.grid(mtry = c(3, 4, 5, 6), 
                        ntree = c(500, 1000), 
                        nodesize = c(30, 50, 70),
                        maxnodes = c(3, 5, 8, 10, 12)
  )
  # split on train/test
  n_repeats <- 30
  
  for (i in seq(1, n_repeats)){
    
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
    data_tr <- limit_outliers(x_train = x_train, x_test = x_test,
                              features = all_features_cols, iqr_multiplicator = 3)
    x_train <- data_tr[['train']]
    x_test <- data_tr[['test']]
    
    # fit models
    trCtrl <- trainControl(
      method="none",
      verboseIter = TRUE,
      classProbs=TRUE,
      summaryFunction = twoClassSummary,
      seeds = seq(1, nrow(par_df)))
    # set "seeds" parameters to make results reproducible
    
    res_iter <- foreach(j=seq(1, nrow(par_df)), .combine = rbind) %dopar% {
      
      mtry <- par_df[j, ]$mtry
      ntree <- par_df[j, ]$ntree
      nodesize <- par_df[j, ]$nodesize
      maxnodes <- par_df[j, ]$maxnodes
      
      mtryGrid <- expand.grid(mtry = mtry)
      model <- train(
        x = x_train, 
        y = y_train,
        preProcess = c("zv", "YeoJohnson", "center", "scale"),
        method = "rf",
        metric="ROC",   
        maximize = TRUE,
        trControl = trCtrl,
        ntree = ntree,
        nodesize = nodesize,
        maxnodes = maxnodes,
        sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
        tuneGrid = mtryGrid
      )
      
      train_pred <- predict(model, newdata = x_train, type = "prob")
      test_pred <- predict(model, newdata = x_test, type = "prob")
      
      tr_auc <- auc(response = y_train, predictor = train_pred[, 2], levels=c("X0", "X1"),  direction = "<") 
      te_auc <- auc(response = y_test, predictor = test_pred[, 2], levels=c("X0", "X1"),  direction = "<") 
      
      out_df <- data.frame(
        mtry = mtry, ntree = ntree, nodesize = nodesize, maxnodes = maxnodes,
        tr_auc = tr_auc, te_auc = te_auc, iter = j)
      
      return(out_df)
    }
    
    res_iter$iter <- i
    res_iter$cancer_type <- cancer_type
    all_results <- rbind(all_results, res_iter)
    
  }
  pb$tick()
  write.csv(all_results, file = "data/output/model_selection/hp_rf_lo.csv")
}

## select not to limit outlies





######## REMOVE CORRELATED FEATURES WITH PLS AND AGAIN TUNE RANDOM FOREST

# https://www.kaggle.com/sasali/notebook83cb3b1069#PLS---RF(customized-caret-train-function)

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

set.seed(7)
tr_neg <- createDataPartition(all_neg_label$strata, p = init_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

all_results <- data.frame()

# set progress bar
n_iterations <- length(all_target_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- data[data[target_column] == 1, ]
  set.seed(7)
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  all_features_cols <- c(conserved_features, grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  par_df <- expand.grid(mtry = c(2, 3, 4), 
                        ntree = c(500, 1000), 
                        nodesize = c(50, 70),
                        maxnodes = c(4, 6, 8),
                        ncomp = c(5, 10, 15)
  )
  
  # split on train/test
  n_repeats <- 30
  
  for (i in seq(1, n_repeats)){
    
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
    # data_tr <- limit_outliers(x_train = x_train, x_test = x_test, 
    #                           features = all_features_cols, iqr_multiplicator = 3)
    # x_train <- data_tr[['train']]
    # x_test <- data_tr[['test']]
    
    # fit models
    trCtrl <- trainControl(
      method="none",
      verboseIter = TRUE,
      classProbs=TRUE,
      summaryFunction = twoClassSummary,
      seeds = seq(1, nrow(par_df)))
    
    res_iter <- foreach(j=seq(1, nrow(par_df)), .combine = rbind) %dopar% {
      
      mtry <- par_df[j, ]$mtry
      ntree <- par_df[j, ]$ntree
      nodesize <- par_df[j, ]$nodesize
      maxnodes <- par_df[j, ]$maxnodes
      ncomp <- par_df[j, ]$ncomp
      
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
        preProcess = c("zv", "YeoJohnson", "center", "scale"),
        method = PlsRf,
        metric="ROC",   
        maximize = TRUE,
        trControl = trCtrl,
        sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
        tuneGrid = mtryGrid
      )
      
      train_pred <- predict(model, newdata = x_train, type = "prob")
      test_pred <- predict(model, newdata = x_test, type = "prob")
      
      tr_auc <- auc(response = y_train, predictor = train_pred[, 2], levels=c("X0", "X1"),  direction = "<") 
      te_auc <- auc(response = y_test, predictor = test_pred[, 2], levels=c("X0", "X1"),  direction = "<")
      
      out_df <- data.frame(
        mtry = mtry, ntree = ntree, nodesize = nodesize, maxnodes = maxnodes, 
        ncomp = ncomp, tr_auc = tr_auc, te_auc = te_auc, iter = j)
      
      return(out_df)
    }
    
    res_iter$iter <- i
    res_iter$cancer_type <- cancer_type
    all_results <- rbind(all_results, res_iter)
    
  }
  pb$tick()
  write.csv(all_results, file = "data/output/model_selection/hp_rf_pls_lo.csv")
}

all_results$diff_auc <- all_results$tr_auc - all_results$te_auc









######## REMOVE CORRELATED FEATURES WITH PCA AND AGAIN TUNE RANDOM FOREST

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

set.seed(7)
tr_neg <- createDataPartition(all_neg_label$strata, p = init_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

all_results <- data.frame()

# set progress bar
n_iterations <- length(all_target_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)


for (target_column in all_target_cols){
  
  set.seed(7)
  cancer_pos_label_data <- data[data[target_column] == 1, ]
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  all_features_cols <- c(conserved_features, grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  par_df <- expand.grid(mtry = c(2, 3, 4), 
                        ntree = c(500, 1000), 
                        nodesize = c(50, 70),
                        maxnodes = c(8, 10, 12),
                        pcaComp = c(5, 10, 15)
  )
  
  # split on train/test
  n_repeats <- 30
  
  for (i in seq(1, n_repeats)){
    
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
    data_tr <- limit_outliers(x_train = x_train, x_test = x_test, 
                              features = all_features_cols, iqr_multiplicator = 3)
    x_train <- data_tr[['train']]
    x_test <- data_tr[['test']]
    
    # fit models
    res_iter <- foreach(j=seq(1, nrow(par_df)), .combine = rbind) %dopar% {
      
      mtry <- par_df[j, ]$mtry
      ntree <- par_df[j, ]$ntree
      nodesize <- par_df[j, ]$nodesize
      maxnodes <- par_df[j, ]$maxnodes
      pcaComp <- par_df[j, ]$pcaComp
      
      mtryGrid <- expand.grid(mtry = mtry)
      model <- train(
        x = x_train, 
        y = y_train,
        preProcess = c("zv", "YeoJohnson", "center", "scale", "pca"),
        method = "rf",
        metric="ROC",   
        maximize = TRUE,
        trControl = trainControl(
          method="none",
          verboseIter = TRUE,
          classProbs=TRUE,
          summaryFunction = twoClassSummary,
          preProcOptions = list(pcaComp = pcaComp),
          seeds = seq(1, nrow(par_df))
        ),
        ntree = ntree,
        nodesize = nodesize,
        maxnodes = maxnodes,
        sampsize = c(X0 = floor(length(tr_pos) * 5), X1 = length(tr_pos)),
        tuneGrid = mtryGrid
      )
      
      train_pred <- predict(model, newdata = x_train, type = "prob")
      test_pred <- predict(model, newdata = x_test, type = "prob")
      
      tr_auc <- auc(response = y_train, predictor = train_pred[, 2], levels=c("X0", "X1"),  direction = "<") 
      te_auc <- auc(response = y_test, predictor = test_pred[, 2], levels=c("X0", "X1"),  direction = "<")
      
      out_df <- data.frame(
        mtry = mtry, ntree = ntree, nodesize = nodesize, maxnodes = maxnodes, pcaComp = pcaComp,
        tr_auc = tr_auc, te_auc = te_auc, iter = j)
      
      return(out_df)
    }
    
    res_iter$iter <- i
    res_iter$cancer_type <- cancer_type
    all_results <- rbind(all_results, res_iter)
    
  }
  pb$tick()
  write.csv(all_results, file = "data/output/model_selection/hp_rf_pca_lo.csv")
}








########################### ADABOOST LO

# discrete adaboost with classification trees, exponential loss 

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

set.seed(7)
tr_neg <- createDataPartition(all_neg_label$strata, p = init_ratio, list = FALSE)
neg_data_train <- all_neg_label[tr_neg, ]

all_results <- data.frame()

for (target_column in all_target_cols){
  
  cancer_pos_label_data <- data[data[target_column] == 1, ]
  set.seed(7)
  tr_pos <- sample(rownames(cancer_pos_label_data), size = floor(init_ratio * nrow(cancer_pos_label_data)), 
                   replace = FALSE)
  pos_data_train <- cancer_pos_label_data[tr_pos, ]
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  all_features_cols <- c(conserved_features, grep(x =  tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  
  par_df <- expand.grid(nu = c(1, 0.1, 2), 
                        iter = c(10, 50, 100, 500), 
                        maxdepth = c(2, 3, 4)
  )  
  
  # split on train/test
  n_repeats <- 30

  for (i in seq(1, n_repeats)){
    
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
    data_tr <- limit_outliers(x_train = x_train, x_test = x_test, 
                              features = all_features_cols, iqr_multiplicator = 3)
    x_train <- data_tr[['train']]
    x_test <- data_tr[['test']]
    
    # fit models
    trCtrl <- trainControl(
      method="none",
      verboseIter = TRUE,
      classProbs=TRUE,
      summaryFunction = twoClassSummary,
      seeds = seq(1, nrow(par_df)))
    
    res_iter <- foreach(j=seq(1, nrow(par_df)), .combine = rbind) %dopar% {

      iter <- par_df[j, ]$iter
      maxdepth <- par_df[j, ]$maxdepth
      nu <- par_df[j, ]$nu
      
      iterGrid <- expand.grid(iter = iter, maxdepth = maxdepth, nu = nu)
      model <- train(
        x = x_train, 
        y = y_train,
        preProcess = c("zv", "YeoJohnson", "center", "scale"),
        method = "ada",
        metric="ROC",   
        maximize = TRUE,
        trControl = trCtrl,
        tuneGrid = iterGrid
      )

      train_pred <- predict(model, newdata = x_train, type = "prob")
      test_pred <- predict(model, newdata = x_test, type = "prob")
      
      tr_auc <- auc(response = y_train, predictor = train_pred[, 2], levels=c("X0", "X1"),  direction = "<") 
      te_auc <- auc(response = y_test, predictor = test_pred[, 2], levels=c("X0", "X1"),  direction = "<")
      
      out_df <- data.frame(
        ntrees = iter, maxdepth = maxdepth, nu = nu, tr_auc = tr_auc, te_auc = te_auc, iter = j)
      return(out_df)
    }
    
    res_iter$iter <- i
    res_iter$cancer_type <- cancer_type
    all_results <- rbind(all_results, res_iter)
    
  }
  
}
all_results$diff_auc  <- NULL
all_results$diff_auc <- all_results$tr_auc - all_results$te_auc
write.csv(all_results, file = "data/output/model_selection/hp_ada_lo.csv")

## analysis of results

all_results <- read.csv("data/output/model_selection/hp_ada_lo.csv")
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

pred_cancers_data <- all_results

ggplot(all_results, aes(te_auc)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~ cancer_type)


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
    ntrees == 100,
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


