script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)

output_folder <- "data/datasets/"



##### Get breakpoints hotspots (target)

# read breakpoints data splitted on 10kb windows
breakpoints_path <- "data/target/"
bkpt <- read.csv(paste0(breakpoints_path, "all_data_10000.csv"), stringsAsFactors = FALSE)
select_cols <- c("chr", "from", "to",  
                 "bkpt_in_window_blood", "bkpt_in_window_bone", "bkpt_in_window_brain",
                 "bkpt_in_window_breast", "bkpt_in_window_liver", "bkpt_in_window_ovary",
                 "bkpt_in_window_pancreatic", "bkpt_in_window_prostate", 
                 "bkpt_in_window_skin", "bkpt_in_window_uterus")
bkpt <- bkpt[select_cols]




win_len <- 100000

# aggregate to win_len
if (win_len == 10000){
  bkpt_data <- bkpt
} else {
  bkpt_data <- bkpt
  bkpt_data$new_from <- floor(bkpt_data$from / win_len) * win_len
  bkpt_data$new_to <- bkpt_data$new_from + win_len
  
  bkpt_cols <- grep(x = names(bkpt_data), pattern = "bkpt", value = TRUE)
  
  
  bkpt_data <- bkpt_data %>%
    group_by(chr, new_from, new_to) %>%
    summarise_at(bkpt_cols, sum)
  
  bkpt_data <- bkpt_data %>%
    rename("to" = "new_to") %>%
    rename("from" = "new_from")
}




# get intersection with good bins
bins_path <- "data/good_bins/"
bins <- read.csv(paste0(bins_path, "bins_", format(win_len, scientific = FALSE), ".csv"), stringsAsFactors = FALSE, 
                 row.names = 1)
chrs <-
  c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )

bins <- bins %>%
  select(chr, start, end) %>%
  setNames(c("chr", "from", "to")) %>% 
  filter(chr %in% chrs)

bkpt_data <- bins %>% 
  inner_join(bkpt_data, by=c("chr", "from", "to"))


# density
quantity_names <- c("bkpt_in_window_blood", "bkpt_in_window_bone", "bkpt_in_window_brain",
                    "bkpt_in_window_breast", "bkpt_in_window_liver", "bkpt_in_window_ovary",
                    "bkpt_in_window_pancreatic", "bkpt_in_window_prostate", 
                    "bkpt_in_window_skin", "bkpt_in_window_uterus")

total_counts <- bkpt_data %>%
  group_by(chr) %>%
  summarise_at(quantity_names, sum) %>%
  ungroup() %>%
  rename_at(quantity_names, function(x) paste0("total_", x))

bkpt_data <- bkpt_data %>%
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
write.csv(densities_id, "data/target/density_100000.csv")

cor_densities <- cor(densities, method="spearman")

for (i in 1:10){
  cor_densities[i, i] <- NA
}
mean(cor_densities, na.rm=TRUE)
median(cor_densities, na.rm=TRUE)
max(cor_densities, na.rm=TRUE)
min(cor_densities, na.rm=TRUE)

# create hotspots
q <- c(0.99, 0.995, 0.999, 0.9995, 0.9999)
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

# check the same target columns and exclude
dupl_cols <- duplicated(t(bkpt_data))
dupl_cols <- names(bkpt_data)[dupl_cols]
duplicates <- list()
for (dupl_col in dupl_cols){
  for (col in setdiff(hsp_cols, c(dupl_cols))){
    ifelse(all(bkpt_data[col] == bkpt_data[dupl_col]), duplicates[dupl_col] <- col, FALSE)
  }
}
duplicates
# win_len <- 10000
# drop_cols <- c("hsp_99%_blood", 
#                "hsp_99%_bone", "hsp_99.5%_bone", 
#                "hsp_99%_brain", "hsp_99.5%_brain")

# win_len <- 20000
# drop_cols <- c("hsp_99%_blood", 
#                "hsp_99%_bone", 
#                "hsp_99%_brain")

# win_len <- 1000000
drop_cols <- c("hsp_99.9%_blood")

bkpt_data <- bkpt_data %>%
  select(-drop_cols)






##### Join features


features_path <- "data/features/"
all_features_paths <- list.files(features_path)


# collect conserved features
conserved_features_path <- c(
  "a-phased_repeats", "direct_repeats", "inverted_repeats",
  "mirror_repeats","short_tandem_repeats",
  "G_quadruplexes",
  "sl_long", "sl_short",
  "genome_regions", "tad","z-dna_motifs"
  )

all_full_conserved_features_path <- vector()

for (folder in conserved_features_path){
  full_path <- paste0(features_path, folder, "/")
  fl <- list.files(full_path)
  all_full_conserved_features_path <- append(
    all_full_conserved_features_path, 
    paste0(
      full_path, 
      grep(x = fl, pattern = paste0("_", format(win_len, scientific = FALSE), ".csv"), 
           value = TRUE)
    )
        )
}

all_data <- bkpt_data

for (feat_path in all_full_conserved_features_path){
  
  feature_data <- read.csv(feat_path, stringsAsFactors = FALSE, header = TRUE)
  if ("X" %in% names(feature_data)){
    feature_data$X <- NULL
  }
  all_data <- all_data  %>%
    left_join(feature_data, by = c("chr", "from", "to"))
  
  feat_col <- names(feature_data)[ncol(feature_data)]
  all_data[is.na(all_data[feat_col]), feat_col] <- 0  
}


# collect tisssue-specific features

tissue_spec_features_path <- c("chromatin_state", "DNA_methylation", "histones", "TF")

all_full_spec_features_path <- vector()

for (folder in tissue_spec_features_path){
  full_path <- paste0(features_path, folder, "/")
  fl <- list.files(full_path)
  all_full_spec_features_path <- append(
    all_full_spec_features_path, 
    paste0(
      full_path, 
      grep(x = fl, pattern = paste0("_", format(win_len, scientific = FALSE), ".csv"), 
           value = TRUE)
    )
  )
}


for (feat_path in all_full_spec_features_path){
  
  feature_data <- read.csv(feat_path, stringsAsFactors = FALSE, header = TRUE)
  if ("X" %in% names(feature_data)){
    feature_data$X <- NULL
  }
  all_cancer_types <- unique(feature_data$cancer_type)
  
  for (canc_type in all_cancer_types){
    
    feature_data_part <- feature_data[feature_data$cancer_type == canc_type, ]
    feature_data_part$cancer_type <- NULL
    feat_name <- names(feature_data_part)[ncol(feature_data_part)]
    new_feat_name <- paste0("cancer_", canc_type, "_", feat_name)
    feature_data_part <- feature_data_part %>%
      rename(!!new_feat_name := feat_name)
    
    all_data <- all_data  %>%
      left_join(feature_data_part, by=c("chr", "from", "to"))
    all_data[is.na(all_data[new_feat_name]), new_feat_name] <- 0  
    
  }
}

write.csv(all_data, 
          file = paste0(output_folder, "dataset_", 
                        format(win_len, scientific = FALSE), ".csv"),
          row.names = FALSE
          )

