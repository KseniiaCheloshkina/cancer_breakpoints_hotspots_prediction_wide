# Boruta considered both highly correlated variables to be important.
# It implies it does not treat collinearity while selecting important variables. It is because of the way algorithm works.

# https://www.listendata.com/2017/05/feature-selection-boruta-package.html
# https://www.analyticsvidhya.com/blog/2016/03/select-important-variables-boruta-package/

# install.packages("Boruta")
# library(Boruta)
# Boruta.Ozone <- Boruta(V4 ~ ., data = Ozone, doTrace = 2, ntree = 1000, mtry = 5)
# не подходит, так как можно регулировать только два параметра random forest, а нам нужно регулировать все

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

# load data
data_path <- "data/datasets/" 
data <- read.csv(
  paste0(data_path, "dataset_", format(win_len, scientific = FALSE), ".csv")
)

output_path <- "data/output/" 

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
hsp_cols <- hsp_cols[-grep(x = hsp_cols, pattern = "all")]

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

n_boruta_resamples <- 10
n_boruta_iter <- 5

# for storing importance results
df_imp_all <- data.frame()

# set progress bar
n_iterations <- nrow(best_lev)
pb <- progress_bar$new(
  format = "  Modeling [:bar] :percent. Elapsed: :elapsedfull ETA: :eta",
  total = n_iterations, clear = FALSE, width=120)
pb$tick(0)




for (target_column in hsp_cols){
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  agg_level <- strsplit(target_column, "_")[[1]][2]
  
  is_best <- best_lev %>%
    inner_join(
      data.frame("cancer_type" = cancer_type, "agg_level" = agg_level),
      by=c("cancer_type","agg_level")
    ) %>% nrow()
  
  if (is_best == 1){
    all_features_cols <- c(all_conserved_feats, grep(x = all_tissue_spec_feats, pattern = cancer_type, value = TRUE))
    
    repeats_res <- foreach(i=seq(1, n_repeats), .combine='rbind') %dopar% {
      
      # train/test split
      splitted_dataset <- get_train_test_split(data=all_data, target_col=target_column, start_pos_col="from",
                                               chr_col="chr", feature_cols=all_features_cols, 
                                               train_ratio=train_ratio, seed=i)
      x_train <- splitted_dataset[["x_train"]]
      y_train <- splitted_dataset[["y_train"]]
      x_test <- splitted_dataset[["x_test"]] 
      y_test <- splitted_dataset[["y_test"]]  
      
      n_pos <- length(y_train[y_train == "X1"])
      n_neg <- length(y_train[y_train == "X0"])
      
      # preprocessing
      preProc  <- preProcess(x_train, method = c("zv", "YeoJohnson", "center", "scale"))
      x_train <- predict(preProc, x_train)
      
      imp_all <- data.frame()
      
      # boruta loop
      for (k in 1:n_boruta_resamples){
        
        iter_features <- names(x_train)
        
        for (j in 1:n_boruta_iter){
          
          # create shadow features
          set.seed(j * k)
          sampl_ind <- sample(x = rownames(x_train), size = length(rownames(x_train)), replace = F)
          x_train_copy <- x_train[sampl_ind, iter_features]
          shadow_names <- paste0("shadow_", names(x_train_copy))
          names(x_train_copy) <- shadow_names
          x_train_new <- cbind(x_train[iter_features], x_train_copy)
          
          # fit model
          rf <- randomForest(x = x_train_new, y = y_train,
                             mtry = 5, ntree = 500, nodesize = 30, maxnodes = 3, 
                             sampsize = get_sample_size(n_pos, n_neg),
                             importance = T)
          
          # take importance
          imp_means <- data.frame(rf$importance)
          imp_means['feature'] <- rownames(imp_means)
          imp_sd <- data.frame(rf$importanceSD)
          imp_sd <- imp_sd %>%
            mutate(
              feature = rownames(imp_sd),
              sd = MeanDecreaseAccuracy
            )
          imp <- imp_means %>%
            select(feature, MeanDecreaseAccuracy) %>%
            inner_join(
              imp_sd %>% select(feature, sd), by = c("feature")
            )
          imp$z_score <- imp$MeanDecreaseAccuracy/imp$sd
          imp[is.na(imp)] <- 0
          
          # select max imp of shadow features
          mzsa <- imp %>%
            filter(feature %in% shadow_names) %>%
            summarize(max(z_score)) %>%
            as.numeric()
          
          res <- imp %>% 
            select(feature, z_score) %>%
            mutate(hit = ifelse(z_score > mzsa, 1, 0),
                   resample = k, iter = j)
          
          iter_features <- res %>%
            filter(!feature %in% shadow_names) %>%
            filter(hit == 1) %>%
            select(feature) 
          
          imp_all <- rbind(imp_all, res)
          
          if (nrow(iter_features) < 5){
            break
          } else {
            iter_features <- iter_features$feature
          }
        }
      }
      imp_all['repeat'] <- i
      
      return(imp_all)
    }
    
    repeats_res <- repeats_res %>%
      mutate(
        cancer_type = cancer_type,
        agg_level = agg_level,
        win_len = format(win_len, scientific = F)
      )
    
    repeats_res$feature_type <- "real"
    repeats_res$feature_type[grep(x = repeats_res$feature, pattern = "shadow")] <- "shadow"
    
    df_imp_all <- rbind(df_imp_all, repeats_res)
    write.csv(df_imp_all, file = paste0(output_path, "boruta_res_", format(win_len, scientific = F), ".csv"))
    
    pb$tick()
    
  }
  
}



