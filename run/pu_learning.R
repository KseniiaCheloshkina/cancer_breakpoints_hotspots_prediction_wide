script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(reshape2)
library(openxlsx)

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

# load data
data_path <- "data/datasets/" 
data <- read.csv(
  paste0(data_path, "dataset_", format(win_len, scientific = FALSE), ".csv")
)

output_path <- "data/output/pu_learning/"


## load best features data
best_features <- read.xlsx("data/output/boruta_plus_sign_stats.xlsx", sheet = "hotspots_features")


# load densities
df_density <- read.csv("data/target/density_100000.csv", stringsAsFactors = F, row.names = 1)

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


# only best labeling type for hotspots
best_lev <- read.csv("data/output/best_lev.csv", row.names = 1)
best_lev$target_name <- paste0("hsp_", best_lev$agg_level, "_", best_lev$cancer_type)
hsp_cols <- best_lev$target_name


# train/test split
train_ratio <- 0.7
n_repeats <- 30

# remove unused columns
all_features_available <- unique(best_features$feature)
all_f <- c(all_conserved_feats, all_tissue_spec_feats)
drop_f <- setdiff(all_f, all_features_available)
all_data <- all_data %>%
  select(-drop_f)

# join density values
all_data <- all_data %>%
  inner_join(
    df_density, by = c("chr", "from", "to")
  )
dens_cols <- grep(x = names(all_data), pattern = "density", value = TRUE)

# for storing results
all_res <- data.frame()
df_recall <- data.frame()

# select quantiles for quality assessment
recall_quantiles <- c(0.001, 0.005, 0.002,  0.003, 0.015, 0.025, seq(0.01, 0.05, 0.01), seq(0.1, 0.9, 0.05))

# set progress bar
n_iterations <- length(hsp_cols)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)

###################### RELIABLE NEGATIVES AND POSITIVES  


for (target_column in hsp_cols){
  
  cat(target_column, file = "data/output/test.log", append = TRUE)
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  agg_level <- strsplit(target_column, "_")[[1]][2]
  dens_col <- grep(x=dens_cols, cancer_type, value=T)
  
  all_features_cols <- c(all_conserved_feats, grep(x = all_tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  # select best features
  features_nm <- best_features[best_features$cancer_type == cancer_type, "feature"]
  
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
  if (length(features_nm) < 5){
    mtryGrid <- expand.grid(mtry = length(features_nm))
  } else {
    mtryGrid <- expand.grid(mtry = 5)
  }
  
  repeats_res <- foreach(i=seq(1, 30)) %dopar% {
    cat("\n", file = "data/output/test.log", append = TRUE)
    cat(paste0("i : ", i), file = "data/output/test.log", append = TRUE)
    
    # train/test split
    splitted_dataset <- get_train_test_split(data=all_data, target_col=target_column, start_pos_col="from",
                                             chr_col="chr", feature_cols=features_nm, train_ratio=train_ratio, 
                                             seed=i, dens_col=dens_col)
    x_train <- splitted_dataset[["x_train"]]
    y_train <- splitted_dataset[["y_train"]]
    x_test <- splitted_dataset[["x_test"]] 
    y_test <- splitted_dataset[["y_test"]]  
    density_test <- splitted_dataset[["density_test"]]  
    
    pos_ind <- which(y_train == "X1")
    n_pos <- length(y_train[y_train == "X1"])
    n_neg <- length(y_train[y_train == "X0"])
  
    epsilon_variants <- c(0.01, 0.03, 0.05)
    n_iter_pu <- 5
    all_test_pred_eps <- data.frame()
    df_recall_all <- data.frame()
    old_rn_i <- 0
    
    cat(paste0("i: ",i), file = "data/output/test.log", append = TRUE)
    cat(paste0("x_train shape: ", nrow(x_train)), file = "data/output/test.log", append = TRUE)
    cat(paste0("x_test shape: ", nrow(x_test)), file = "data/output/test.log", append = TRUE)
    df_agg_stats_eps <- data.frame()
    
    # perform pu learning
    for (eps in epsilon_variants){
      print(eps)
      cat(paste0("eps : ", eps), file = "data/output/test.log", append = TRUE)
      # cat(paste0("eps: ", eps), file = "data/output/test.log", append = TRUE)
      
      # initial step
      all_test_pred <- data.frame()
      new_x_train <- x_train
      new_y_train <- y_train
      
      # fit model
      print('nPOS')
      print(n_pos)
      print('nNEG')
      print(n_neg)
      cat(paste0("n_pos: ",n_pos), file = "data/output/test.log", append = TRUE)
      cat(paste0("n_neg: ",n_neg ), file = "data/output/test.log", append = TRUE)
      
      model <- rf_fit(x_train = new_x_train, y_train = new_y_train, trCtrl = trCtrl, 
                      n_pos = n_pos, n_neg = n_neg, 
                      mtryGrid = mtryGrid)
      
      # get predictions
      train_pred <- predict(model, newdata = new_x_train, type = "prob")
      train_pred$target <- new_y_train
      
      test_pred <- predict(model, newdata = x_test, type = "prob")
      test_pred$target <- y_test
      
      # model quality
      model_qual <- get_model_quality(train_pred, test_pred, model, recall_quantiles)
      model_qual <- model_qual[['recall']]
      model_qual$eps <- eps
      model_qual$iter_pu <- 0
      
      cat(paste0("new_x_train shape: ", nrow(new_x_train)), file = "data/output/test.log", append = TRUE)
      cat(paste0("x_test shape: ", nrow(x_test)), file = "data/output/test.log", append = TRUE)

      # add density to predictions
      test_pred$density <-  density_test[, 1]
      
      # cutoff
      q1 <- quantile(train_pred[train_pred$target == "X1", "X1"], 0.1)
      q9 <- quantile(train_pred[train_pred$target == "X1", "X1"], 0.9)
      if (q1 == q9) {
        next
      } else if (q9 - q1 <= 2 * eps){
        rlow <- q1
        rup <- q9
      } else {
        rlow <- q1 + eps
        rup <- q9 - eps
      }
      cat(paste0("rlow : ", rlow), file = "data/output/test.log", append = TRUE)
      cat(paste0("rup : ", rup), file = "data/output/test.log", append = TRUE)
      
      # test labeling
      test_pred$rp <- 0
      test_pred[test_pred$X1 >= rup, "rp"] <- 1
      test_pred$iter <- 0
      test_pred$rn <- 0

      all_test_pred <- rbind(all_test_pred, test_pred)

      te_auc <- as.numeric(
        auc(cases = test_pred$X1[test_pred$target == "X1"],  controls = test_pred$X1[test_pred$target == "X0"],
            direction = "<")
      )
      
      df_agg_stats <- data.frame(
        iter=0,
        train_total=nrow(x_train),
        train_rn=length(which(y_train == "X0")),
        train_rp=length(which(y_train == "X1")),
        train_unl=0,
        roc_auc_test=te_auc
      )      

      # reliable negatives and positives
      all_ind <- which(train_pred$target == "X0")
      rn_ind <- which(train_pred[train_pred$target == "X0", "X1"] < rlow)
      rp_ind <- which(train_pred[train_pred$target == "X0", "X1"] > rup)
      unl_ind <- setdiff(all_ind, c(rn_ind, rp_ind))
      
      rn_i <- length(rn_ind)
      rp_i <- length(rp_ind)
      q_i <- length(unl_ind)
      
      cat(paste0("rn_i : ", rn_i), file = "data/output/test.log", append = TRUE)
      cat(paste0("q_i : ", q_i), file = "data/output/test.log", append = TRUE)
      
      all_test_pred$eps <- eps
      all_test_pred_eps <- rbind(all_test_pred_eps, all_test_pred)
      df_agg_stats$eps <- eps
      df_agg_stats_eps <- rbind(df_agg_stats_eps, df_agg_stats)
      df_recall_all <- rbind(df_recall_all, model_qual)
      
      if ((q_i == 0) | (rn_i == 0)){
        next
      } 
      
      # new labeled and unlabeled dataset
      unl_x_train <- new_x_train[unl_ind, ]
      unl_y_train <- new_y_train[unl_ind]
      
      new_x_train <- new_x_train[c(rp_ind, pos_ind, rn_ind), ]
      new_y_train <- as.factor(c(rep("X1", rp_i), as.character(new_y_train[c(pos_ind, rn_ind)])))
      
      cat(paste0("unl_x_train : ", nrow(unl_x_train)), file = "data/output/test.log", append = TRUE)
      cat(paste0("new_x_train : ", nrow(new_x_train)), file = "data/output/test.log", append = TRUE)
  
      for (n in 1:n_iter_pu){
        cat(paste0("n : ", n), file = "data/output/test.log", append = TRUE)
        # fit model
        print('n')
        print(n)
        rn_ind_c <- length(which(new_y_train == "X0"))
        cat(paste0("n_pos: ", length(pos_ind)), file = "data/output/test.log", append = TRUE)
        cat(paste0("n_neg: ", rn_ind_c ), file = "data/output/test.log", append = TRUE)
        
        model <- rf_fit(new_x_train, new_y_train, trCtrl, 
                        n_pos = length(pos_ind), n_neg = rn_ind_c, mtryGrid)
        
        # get predictions
        train_pred <- predict(model, newdata = new_x_train, type = "prob")
        train_pred$target <- new_y_train
        test_pred <- predict(model, newdata = x_test, type = "prob")
        test_pred$target <- y_test
        unl_pred <- predict(model, newdata = unl_x_train, type = "prob")

        # model quality
        model_qual <- get_model_quality(train_pred, test_pred, model, recall_quantiles)
        model_qual <- model_qual[['recall']]
        model_qual$eps <- eps
        model_qual$iter_pu <- n
        df_recall_all <- rbind(df_recall_all, model_qual)
        
        # add density to predictions
        test_pred$density <- density_test[, 1]
        
        # cutoff
        q1 <- quantile(train_pred[train_pred$target == "X1", "X1"], 0.1)
        q9 <- quantile(train_pred[train_pred$target == "X1", "X1"], 0.9)

        if (q1 == q9) {
          break
        } else if (q9 - q1 <= 2 * eps){
          rlow <- q1
          rup <- q9
        } else {
          rlow <- q1 + eps
          rup <- q9 - eps
        }        
        
        cat(paste0("rlow : ", rlow), file = "data/output/test.log", append = TRUE)
        cat(paste0("rup : ", rup), file = "data/output/test.log", append = TRUE)
        
        # test labeling
        test_pred$rp <- 0
        test_pred[test_pred$X1 >= rup, "rp"] <- 1
        test_pred$iter <- n
        
        rownames(all_test_pred) <- rownames(test_pred) <- NULL
        test_pred$rn <- rn_i
        test_pred$eps <- eps
        all_test_pred_eps <- rbind(all_test_pred_eps, test_pred)
        
        te_auc <- as.numeric(
          auc(cases = test_pred$X1[test_pred$target == "X1"],  controls = test_pred$X1[test_pred$target == "X0"], direction = "<")
        )
        
        df_agg_stats_eps <- rbind(df_agg_stats_eps, data.frame(
          iter=n,
          train_total=nrow(new_x_train),
          train_rn=length(which(new_y_train == "X0")),
          train_rp=length(which(new_y_train == "X1")),
          train_unl=nrow(unl_x_train),
          roc_auc_test=te_auc,
          eps=eps
        ))
  
        # reliable negatives from unlabeled
        rn_ind <- which(unl_pred$X1 < rlow)
        rp_ind <- which(unl_pred$X1 > rup)
        unl_ind <- which((unl_pred$X1 >= rlow) & (unl_pred$X1 <= rup))
        rn_i <- length(rn_ind)
        rp_i <- length(rp_ind)
        q_i <- length(unl_ind)
        
        cat(paste0("rn_i : ", rn_i), file = "data/output/test.log", append = TRUE)
        cat(paste0("q_i : ", q_i), file = "data/output/test.log", append = TRUE)
        
        # condition for leaving
        # number of labeled examples on one iteration should decline with time
        # number of labeled examples on one iteration should be significant (more than there are positive examples)
        print('conditions for leave')
        print(n)
        print(rn_i)
        print(old_rn_i)
        print(n_pos)
        if ((n > 1) & ((rn_i > old_rn_i) | (rn_i < n_pos))){
          break
        } else {
          
          if ((rup <= rlow) | (q_i == 0) | (rp_i + rn_i == 0)) {
            print(rup)
            print(rlow)
            print(q_i)
            break
          } else {
          
            old_rn_i <- rn_i
            
            # new labeled dataset
            new_x_train <- rbind(new_x_train, unl_x_train[c(rn_ind, rp_ind), ])
            
            new_y_train <- c(as.character(new_y_train), rep("X0", rn_i), rep("X1", rp_i))
            new_y_train <- as.factor(new_y_train)

            unl_x_train <- unl_x_train[unl_ind, ]
          }
        }

      }

    }
    
    res <- all_test_pred_eps %>% 
      group_by(eps, iter, target) %>% 
      summarize(
        rp = sum(rp)
      ) %>%
      dcast(eps + iter ~ target, value.var = "rp") %>%
      setNames(c("eps", "iter", "rp_0", "rp_hsp")) %>%
      mutate(
        rp_total = rp_0 + rp_hsp
      ) %>%
      inner_join(
        all_test_pred_eps %>% 
          filter(rp == 1) %>%
          group_by(eps, iter) %>% 
          summarize(
            mean_bkpt_density = mean(density)
          ),
        by=c("eps", "iter")
      ) %>%
      select(-rp_0)  %>%
      inner_join(
        all_test_pred_eps %>% 
          group_by(eps, iter) %>% 
          summarize(
            mean_bkpt_density_test = mean(density)
          ),
        by=c("eps", "iter")
      ) %>%
      inner_join(
        df_agg_stats_eps,
        by = c("eps", 'iter')
      ) %>%
      mutate(
        test_total = length(y_test),
        test_hsp = length(which(y_test == "X1")),
        resample = i,
        all_train_possible = nrow(x_train)
      ) 
    
    final_res <- list()
    final_res[['stats']] <- res
    final_res[['recall']] <- df_recall_all
    
    return(final_res)
  }

  # save results
  for (split_iter in 1:length(repeats_res)){
    
    res_iter <- repeats_res[[split_iter]]
    
    # extract specific datasets
    df_recall_all <- res_iter[['recall']]
    stats <- res_iter[['stats']]  
    
    df_recall_all <- df_recall_all %>%
      mutate(
        iter = split_iter,
        cancer_type = cancer_type,
        agg_level = agg_level,
        win_len = format(win_len, scientific = F)
      )
    df_recall <- rbind(df_recall, df_recall_all)
    
    stats$cancer_type <- cancer_type
    stats$agg_level <- agg_level
    all_res <- rbind(all_res, stats)
  }
  
  write.csv(df_recall, paste0(output_path, "pu_learning_rp_recall.csv"))
  write.csv(all_res, paste0(output_path, "pu_learning_rp.csv"))
  pb$tick()
}


###################### RELIABLE NEGATIVES ONLY

for (target_column in hsp_cols){
  
  cat(target_column, file = "data/output/test.log", append = TRUE)
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  agg_level <- strsplit(target_column, "_")[[1]][2]
  dens_col <- grep(x=dens_cols, cancer_type, value=T)
  
  all_features_cols <- c(all_conserved_feats, grep(x = all_tissue_spec_feats, pattern = cancer_type, value = TRUE))
  
  # select best features
  features_nm <- best_features[best_features$cancer_type == cancer_type, "feature"]
  
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
  if (length(features_nm) < 5){
    mtryGrid <- expand.grid(mtry = length(features_nm))
  } else {
    mtryGrid <- expand.grid(mtry = 5)
  }
  
  repeats_res <- foreach(i=seq(1, n_repeats)) %dopar% {
    cat("\n", file = "data/output/test.log", append = TRUE)
    cat(paste0("i : ", i), file = "data/output/test.log", append = TRUE)
    
    # train/test split
    splitted_dataset <- get_train_test_split(data=all_data, target_col=target_column, start_pos_col="from",
                                             chr_col="chr", feature_cols=features_nm, train_ratio=train_ratio, 
                                             seed=i, dens_col=dens_col)
    x_train <- splitted_dataset[["x_train"]]
    y_train <- splitted_dataset[["y_train"]]
    x_test <- splitted_dataset[["x_test"]] 
    y_test <- splitted_dataset[["y_test"]]  
    density_test <- splitted_dataset[["density_test"]]  
    
    pos_ind <- which(y_train == "X1")
    n_pos <- length(y_train[y_train == "X1"])
    n_neg <- length(y_train[y_train == "X0"])
    
    epsilon_variants <- c(0.01, 0.03, 0.05)
    n_iter_pu <- 5
    all_test_pred_eps <- data.frame()
    df_recall_all <- data.frame()
    old_rn_i <- 0
    
    cat(paste0("i: ",i), file = "data/output/test.log", append = TRUE)
    cat(paste0("x_train shape: ", nrow(x_train)), file = "data/output/test.log", append = TRUE)
    cat(paste0("x_test shape: ", nrow(x_test)), file = "data/output/test.log", append = TRUE)
    df_agg_stats_eps <- data.frame()
    
    # perform pu learning
    for (eps in epsilon_variants){
      
      # cat(paste0("eps: ", eps), file = "data/output/test.log", append = TRUE)
      
      # initial step
      all_test_pred <- data.frame()
      new_x_train <- x_train
      new_y_train <- y_train
      
      # fit model
      print('nPOS')
      print(n_pos)
      print('nNEG')
      print(n_neg)
      cat(paste0("n_pos: ",n_pos), file = "data/output/test.log", append = TRUE)
      cat(paste0("n_neg: ",n_neg ), file = "data/output/test.log", append = TRUE)
      
      model <- rf_fit(x_train = new_x_train, y_train = new_y_train, trCtrl = trCtrl, 
                      n_pos = n_pos, n_neg = n_neg, 
                      mtryGrid = mtryGrid)
      
      # get predictions
      train_pred <- predict(model, newdata = new_x_train, type = "prob")
      train_pred$target <- new_y_train
      
      test_pred <- predict(model, newdata = x_test, type = "prob")
      test_pred$target <- y_test
      
      # model quality
      model_qual <- get_model_quality(train_pred, test_pred, model, 
                                      recall_quantiles)
      model_qual <- model_qual[['recall']]
      model_qual$eps <- eps
      model_qual$iter_pu <- 0
      
      # cat(paste0("new_x_train shape: ", nrow(new_x_train)), file = "data/output/test.log", append = TRUE)
      # cat(paste0("x_test shape: ", nrow(x_test)), file = "data/output/test.log", append = TRUE)
      
      # add density to predictions
      test_pred$density <-  density_test[, 1]
      
      # cutoff
      q1 <- quantile(train_pred[train_pred$target == "X1", "X1"], 0.1)
      q9 <- quantile(train_pred[train_pred$target == "X1", "X1"], 0.9)
      if (q1 == q9) {
        next
      } else if (q9 - q1 <= 2 * eps){
        rlow <- q1
        rup <- q9
      } else {
        rlow <- q1 + eps
        rup <- q9 - eps
      }
      cat(paste0("rlow : ", rlow), file = "data/output/test.log", append = TRUE)
      cat(paste0("rup : ", rup), file = "data/output/test.log", append = TRUE)
      
      # test labeling
      test_pred$rp <- 0
      test_pred[test_pred$X1 >= rup, "rp"] <- 1
      test_pred$iter <- 0
      test_pred$rn <- 0
      all_test_pred <- rbind(all_test_pred, test_pred)
      
      te_auc <- as.numeric(
        auc(cases = test_pred$X1[test_pred$target == "X1"],  controls = test_pred$X1[test_pred$target == "X0"],
            direction = "<")
      )
      
      df_agg_stats <- data.frame(
        iter=0,
        train_total=nrow(x_train),
        train_rn=length(which(y_train == "X0")),
        train_rp=length(which(y_train == "X1")),
        train_unl=0,
        roc_auc_test=te_auc
      )
      
      # reliable negatives
      all_ind <- which(train_pred$target == "X0")
      rn_ind <- which(train_pred[train_pred$target == "X0", "X1"] < rlow)
      unl_ind <- setdiff(all_ind, rn_ind)
      
      rn_i <- length(rn_ind)
      q_i <- length(unl_ind)
      
      # cat(paste0("rn_i : ", rn_i), file = "data/output/test.log", append = TRUE)
      # cat(paste0("q_i : ", q_i), file = "data/output/test.log", append = TRUE)
      
      all_test_pred$eps <- eps
      all_test_pred_eps <- rbind(all_test_pred_eps, all_test_pred)
      df_recall_all <- rbind(df_recall_all, model_qual)
      df_agg_stats$eps <- eps
      df_agg_stats_eps <- rbind(df_agg_stats_eps, df_agg_stats)
      
      if ((q_i == 0) | (rn_i == 0)){
        next
      }     
      
      # new labeled and unlabeled dataset
      unl_x_train <- new_x_train[unl_ind, ]
      unl_y_train <- new_y_train[unl_ind]
      
      new_x_train <- new_x_train[c(pos_ind, rn_ind), ]
      new_y_train <- new_y_train[c(pos_ind, rn_ind)] 
      
      # cat(paste0("unl_x_train : ", nrow(unl_x_train)), file = "data/output/test.log", append = TRUE)
      # cat(paste0("new_x_train : ", nrow(new_x_train)), file = "data/output/test.log", append = TRUE)
      
      for (n in 1:n_iter_pu){

        rn_ind_c <- length(which(new_y_train == "X0"))
        cat(paste0("n_pos: ", length(pos_ind)), file = "data/output/test.log", append = TRUE)
        cat(paste0("n_neg: ", rn_ind_c ), file = "data/output/test.log", append = TRUE)
        # fit model
        model <- rf_fit(new_x_train, new_y_train, trCtrl, 
                        n_pos = n_pos, n_neg = rn_ind_c, mtryGrid)
        
        # get predictions
        train_pred <- predict(model, newdata = new_x_train, type = "prob")
        train_pred$target <- new_y_train
        test_pred <- predict(model, newdata = x_test, type = "prob")
        test_pred$target <- y_test
        unl_pred <- predict(model, newdata = unl_x_train, type = "prob")
        
        # model quality
        model_qual <- get_model_quality(train_pred, test_pred, model, recall_quantiles)
        model_qual <- model_qual[['recall']]
        model_qual$eps <- eps
        model_qual$iter_pu <- n
        df_recall_all <- rbind(df_recall_all, model_qual)
        
        # add density to predictions
        test_pred$density <- density_test[, 1]
        
        # cutoff
        q1 <- quantile(train_pred[train_pred$target == "X1", "X1"], 0.1)
        q9 <- quantile(train_pred[train_pred$target == "X1", "X1"], 0.9)
        
        if (q1 == q9) {
          break
        } else if (q9 - q1 <= 2 * eps){
          rlow <- q1
          rup <- q9
        } else {
          rlow <- q1 + eps
          rup <- q9 - eps
        }        
        
        # cat(paste0("rlow : ", rlow), file = "data/output/test.log", append = TRUE)
        # cat(paste0("rup : ", rup), file = "data/output/test.log", append = TRUE)  
        
        # test labeling
        test_pred$rp <- 0
        test_pred[test_pred$X1 >= rup, "rp"] <- 1
        test_pred$iter <- n
        
        rownames(all_test_pred) <- rownames(test_pred) <- NULL
        test_pred$rn <- rn_i
        test_pred$eps <- eps
        all_test_pred_eps <- rbind(all_test_pred_eps, test_pred)
        
        te_auc <- as.numeric(
          auc(cases = test_pred$X1[test_pred$target == "X1"],  controls = test_pred$X1[test_pred$target == "X0"],
              direction = "<")
        )
        
        df_agg_stats_eps <- rbind(df_agg_stats_eps, data.frame(
          iter=n,
          train_total=nrow(new_x_train),
          train_rn=length(which(new_y_train == "X0")),
          train_rp=length(which(new_y_train == "X1")),
          train_unl=nrow(unl_x_train),
          roc_auc_test=te_auc,
          eps=eps
        ))
        
        # reliable negatives from unlabeled
        rn_ind <- which(unl_pred$X1 < rlow)
        unl_ind <- which(unl_pred$X1 >= rlow)
        
        rn_i <- length(rn_ind)
        q_i <- length(unl_ind)
        # cat(paste0("rn_i : ", rn_i), file = "data/output/test.log", append = TRUE)
        # cat(paste0("q_i : ", q_i), file = "data/output/test.log", append = TRUE)          
        
        # condition for leaving
        # number of labeled examples on one iteration should decline with time
        # number of labeled examples on one iteration should be significant (more than there are positive examples)
        if ((n > 1) & ((rn_i > old_rn_i) | (rn_i < n_pos))){
          break
        } else {
          
          if ((rup <= rlow) | (q_i == 0) | (rn_i == 0)) {
            break
          } else {
            
            old_rn_i <- rn_i
            
            # new labeled dataset
            new_x_train <- rbind(new_x_train, unl_x_train[rn_ind, ])
            
            new_y_train <- c(as.character(new_y_train), as.character(unl_y_train[rn_ind]))
            new_y_train <- as.factor(new_y_train)
            
            unl_x_train <- unl_x_train[unl_ind, ]
          }
        }
      }
    }
    
    res <- all_test_pred_eps %>% 
      group_by(eps, iter, target) %>% 
      summarize(
        rp = sum(rp)
      ) %>%
      dcast(eps + iter ~ target, value.var = "rp") %>%
      setNames(c("eps", "iter", "rp_0", "rp_hsp")) %>%
      mutate(
        rp_total = rp_0 + rp_hsp
      ) %>%
      inner_join(
        all_test_pred_eps %>% 
          filter(rp == 1) %>%
          group_by(eps, iter) %>% 
          summarize(
            mean_bkpt_density = mean(density)
          ),
        by=c("eps", "iter")
      ) %>%
      select(-rp_0)  %>%
      inner_join(
        all_test_pred_eps %>% 
          group_by(eps, iter) %>% 
          summarize(
            mean_bkpt_density_test = mean(density)
          ),
        by=c("eps", "iter")
      ) %>%
      inner_join(
        df_agg_stats_eps,
        by = c("eps", 'iter')
      ) %>%
      mutate(
        test_total = length(y_test),
        test_hsp = length(which(y_test == "X1")),
        resample = i,
        all_train_possible = nrow(x_train)
      ) 
    
    final_res <- list()
    final_res[['stats']] <- res
    final_res[['recall']] <- df_recall_all
    
    return(final_res)
  }
  
  # save results
  for (split_iter in 1:length(repeats_res)){
    
    res_iter <- repeats_res[[split_iter]]
    
    # extract specific datasets
    df_recall_all <- res_iter[['recall']]
    stats <- res_iter[['stats']]  
    
    df_recall_all <- df_recall_all %>%
      mutate(
        iter = split_iter,
        cancer_type = cancer_type,
        agg_level = agg_level,
        win_len = format(win_len, scientific = F)
      )
    df_recall <- rbind(df_recall, df_recall_all)
    
    stats$cancer_type <- cancer_type
    stats$agg_level <- agg_level
    all_res <- rbind(all_res, stats)
  }
  write.csv(df_recall, paste0(output_path, "pu_learning_rn_recall.csv"))
  write.csv(all_res, paste0(output_path, "pu_learning_rn.csv"))
  pb$tick()
}

