library(dplyr)
library(reshape2)
library(stringr)
require(data.table)
library(openxlsx)

script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)

source("../cbp_data/get_intersections.R")

win_len <- 100000

####################### BREAKPOINTS IN LOW-MAPPABLE REGIONS 
# load breakpoints
df <- read.csv("../../cbp_data/data/raw breakpoints/breast_all_data.csv")
df <- df %>% 
  select(chr, chr_bkpt_beg, chr_bkpt_end) %>%
  setNames(c("chr", "start", "end"))

# load blacklist

# The Duke Uniqueness track 35mer
# Scores are normalized to between 0 and 1, with 1 representing a 
# completely unique sequence and 0 representing a sequence that occurs
# more than 4 times in the genome (excluding chrN_random and alternative 
# haplotypes). A score of 0.5 indicates the sequence occurs exactly twice,
# likewise 0.33 for three times and 0.25 for four times.

# select regions with uniqueness score lower than threshold 'thr'
# these regions are called 'low-mappable' (are met more frequently than threshold)
# and should be removed

low_map_path <-'../../cbp_data/data/genome/low_mappability.bedGraph'
df_low_map <- read.table(low_map_path, header=FALSE)
names(df_low_map) <- c("chr", "start", "end", "score")
df_low_map$chr <- as.character(df_low_map$chr)
df_low_map <- df_low_map[!df_low_map$chr %in% c("chrM", "chrY"), ]
df_low_map$chr <- sapply(df_low_map$chr, function(x) gsub(pattern = "chr", replacement = "", x=x))
df_low_map$score <- round(df_low_map$score, 2)

# find breakpoints that lies in low-mappable regions

# table(df_low_map$score)
# 0     0.25 0.333333      0.5        1 
# 9976947  5850063  7685716 11875404 1229957

# cumulative
thrs <- c(0.25, 0.33, 0.5, 1)
all_data <- data.frame()

for (thr in thrs){
  df_low_map_thr <- df_low_map[df_low_map$score < thr, ]
  df_low_map_thr$score <- NULL
  
  windows_with_blacklisted <- get_intersection_intervals(df, df_low_map_thr)
  windows_with_blacklisted <- windows_with_blacklisted %>%
    select(chr, start, end) %>%
    mutate(
      thr = thr
    )
  all_data <- rbind.data.frame(all_data, windows_with_blacklisted)
}
all_data <- all_data %>% unique()
all_data$ones <- 1
all_data_wide <- dcast(all_data, chr + start + end ~ thr, value.var="ones")
all_data_wide[is.na(all_data_wide)] <- 0


stats <- all_data %>%
  unique() %>%
  group_by(thr) %>%
  summarize(
    n_breaks_in_non_unique = n()
  ) %>% 
  mutate(
    proportion_breaks_in_non_unique = n_breaks_in_non_unique/nrow(df)
  )


# select threshold <0.5 and save final breakpoints
df_low_map_thr <- df_low_map[df_low_map$score < 0.5, ]
df_low_map_thr$score <- NULL

windows_with_blacklisted <- get_intersection_intervals(df, df_low_map_thr)
windows_with_blacklisted <- windows_with_blacklisted %>%
  select(chr, start, end) %>%
  mutate(
    intesects_blacklisted = 1
  ) %>% 
  unique()

df <- df %>%
  left_join(windows_with_blacklisted)
df[is.na(df)] <- 0
write.csv(df, file="../data/output_third/reports/breast_low_mappable_data.csv")

# wb <- createWorkbook()
# addWorksheet(wb, "breast_summary")
# addWorksheet(wb, "breast_labels")
# 
# writeData(wb, sheet="breast_summary", stats)
# writeData(wb, sheet="breast_labels", df)
# 
# saveWorkbook(wb, "../data/output_third/reports/breast_low_mappable.xlsx", overwrite = T)







####################### CALCULATE NEW TARGET 


bkpt_data <- read.csv("../data/output_third/reports/breast_low_mappable_data.csv")
bkpt_data <- bkpt_data %>%
  filter(intesects_blacklisted == 0) %>%
  select(chr, start, end) %>%
  mutate(cancer_type = "breast")

# read bins
bins_path <- "../data/good_bins/"
bins <- read.csv(paste0(bins_path, "bins_", format(win_len, scientific = FALSE), ".csv"), 
                 stringsAsFactors = FALSE, row.names = 1)

# get intersection with good bins
all_bkpt_counts <- bins
all_bkpt_counts <- all_bkpt_counts %>%
  select(chr, start, end) %>%
  setNames(c("chr", "from", "to"))  

for (canc in unique(bkpt_data$cancer_type)){
  data_cancer <- bkpt_data[bkpt_data$cancer_type == canc, ]
  data_cancer$cancer_type <- NULL
  data_intersected <- get_intersection_intervals(bins, data_cancer)
  data_intersected <- data_intersected %>%
    group_by(chr, start, end) %>%
    summarize(bkpt_in_window = n())
  
  data_intersected <- data_intersected %>%
    select(chr, start, end, bkpt_in_window) %>%
    setNames(c("chr", "from", "to", paste0("bkpt_in_window_", canc)))
  
  all_bkpt_counts <- all_bkpt_counts %>%
    left_join(data_intersected, by = c('chr', "from", 'to'))
}

all_bkpt_counts[is.na(all_bkpt_counts)] <- 0
all_bkpt_counts <- all_bkpt_counts[all_bkpt_counts$chr != "Y", ]

# density
quantity_names <- c("bkpt_in_window_breast")

total_counts <- all_bkpt_counts %>%
  group_by(chr) %>%
  summarise_at(quantity_names, sum) %>%
  ungroup() %>%
  rename_at(quantity_names, function(x) paste0("total_", x))

bkpt_data <- all_bkpt_counts %>%
  inner_join(total_counts, by=c("chr"))


for (col in quantity_names){
  new_col_name <- paste0("density_", tail(strsplit(x = col, split = "_")[[1]], 1))
  bkpt_data <- bkpt_data %>%
    mutate(!!new_col_name := !!as.name(col) / !!as.name(paste0("total_", col)))
}


densities <- bkpt_data %>%
  select_at(vars(contains("density")))

densities_id <- bkpt_data %>%
  select(c("chr", "from", "to", names(densities)))

# create hotspots
q <- c(0.99)
density_cols <- grep(x = names(bkpt_data), pattern = "density", value = TRUE)

for (col in density_cols){
  q_value <- quantile(x = bkpt_data[, col], probs = q)
  q_vars <- names(q_value)
  for (q_var_cur in q_vars){
    new_col_name <- paste0("hsp_", format(q_var_cur, scientific = FALSE), "_", strsplit(x = col, split = "_")[[1]][2])
    bkpt_data <- bkpt_data %>%
      mutate(!!new_col_name := ifelse(!!as.name(col) > q_value[q_var_cur], 1, 0))    
  }
  # add variant when all windows with non-zero density are labeled as ones
  new_col_name <- paste0("hsp_all_", strsplit(x = col, split = "_")[[1]][2])
  bkpt_data <- bkpt_data %>%
    mutate(!!new_col_name := ifelse(!!as.name(col) > 0, 1, 0))   
}

hsp_cols <- grep(x = names(bkpt_data), pattern = "hsp", value = TRUE)
select_cols <- c("chr", "from", "to", hsp_cols)
bkpt_data <- bkpt_data[select_cols]

# join features
df_dataset <- read.csv("../data/datasets/dataset_100000.csv")
df_dataset$X <- NULL
select_cols <- names(df_dataset)[-grep(names(df_dataset), pattern = "hsp")]
df_dataset <- df_dataset[, select_cols]

df_dataset <- df_dataset %>%
  inner_join(bkpt_data)

write.csv(df_dataset, "../data/datasets/dataset_100000_breast_lm.csv", row.names = FALSE)







####################### BORUTA FEATURE SELECTION
# look at script run/boruta_fs.R

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

n_cores <- 1
registerDoParallel(n_cores)

set.seed(7)


win_len <- 100000
win_len_upper <- win_len * 10


# load data
data_path <- "data/datasets/" 
data <- read.csv("data/datasets/dataset_100000_breast_lm.csv")

output_path <- "data/output_third/" 

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

# hsp_cols <- hsp_cols[-grep(x = hsp_cols, pattern = "all")]

# train/test split
train_ratio <- 0.7
n_repeats <- 30

n_boruta_resamples <- 10
n_boruta_iter <- 5

# for storing importance results
df_imp_all <- data.frame()

for (target_column in hsp_cols){
  
  cancer_type <- strsplit(target_column, "_")[[1]][3]
  agg_level <- strsplit(target_column, "_")[[1]][2]
  
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
  write.csv(df_imp_all, file = paste0(output_path, "boruta_res_", "lm", ".csv"))
  
}








