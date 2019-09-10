script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(tcltk)
library(purrr)
library(foreach)
library(doParallel)

n_cores <- 10
registerDoParallel(n_cores)

source("cbp_data/get_intersections.R")


################## function for calculation of coverage  
################## for overlapped intervals in one window
calculate_coverage_for_overlapped <- function(data, int_start, int_end, win_len){
  data_part <- data[(data$interval_start == int_start) & (data$interval_end == int_end),]
  n_int <- nrow(data_part)
  
  tot_len <- data_part$start[1] - int_start
  start_pos <- data_part$start[1]
  end_pos <- data_part$end[1]
  
  for (i in 2:n_int){
    if (data_part$start[i] < end_pos){
      if (data_part$end[i] > end_pos){
        end_pos <- data_part$end[i]
      }
    }
    else {
      tot_len <- tot_len + (data_part$start[i] - end_pos)
      start_pos <- data_part$start[i]
      end_pos <- data_part$end[i]
    }
  }
  if (int_end > end_pos){
    tot_len <- tot_len + int_end - end_pos
  }
  cov <- (win_len - tot_len) / win_len
  
  return(cov)
}

################## function for calculation of coverage  
################## for not overlapping intervals in all windows
get_coverage_not_overlapped <- function(data_intersected, win_len, feat_name){
  
  data_grouped <- data_intersected %>%
    mutate(length = corr_end - corr_start) %>%
    group_by(chr, from, to) %>%
    summarize(summed_len = sum(length)) %>%
    mutate(cov = summed_len / win_len)
  
  all_data <- data_grouped %>%
    ungroup() %>%
    select(chr, from, to, cov) %>% 
    setNames(c("chr", "from", "to", feat_name))
  
  return(all_data)
}

################## function for calculation of coverage  
################## for overlapping intervals in all windows
get_coverage_overlapped <- function(res, data_intersected, win_len, feat_name){
  
  overlapped_int <- res %>%
    select(chr, from, to) %>%
    unique() %>%
    mutate(
      overlapped = 1
    )
  
  all_int <- data_intersected %>%
    left_join(overlapped_int, by=c("chr", 'from', 'to')) 
  
  overlapped_int <- all_int %>%
    filter(overlapped == 1)
  
  # not overlapped part
  not_overlapped_int <- all_int %>%
    filter(is.na(overlapped))
  
  data_grouped <- not_overlapped_int %>%
    mutate(length = corr_end - corr_start) %>%
    group_by(chr, from, to) %>%
    summarize(summed_len = sum(length)) %>%
    mutate(cov = summed_len / win_len)
  
  data_grouped <- data_grouped %>%
    ungroup() %>%
    select(chr, from, to, cov) %>% 
    setNames(c("chr", "from", "to", feat_name))
  
  # overlapped part
  overlapped_int <- overlapped_int %>%
    select(chr, from, to, corr_start, corr_end) %>% 
    setNames(c("chr", "interval_start", "interval_end",  "start", "end"))
  
  unique_combs <- unique(overlapped_int[, c("chr", "interval_start", "interval_end")])
  
  data_overlapped <- foreach(i=seq(1, nrow(unique_combs)), .combine = rbind) %dopar% {
    cur_chr <- unique_combs[i, ]$chr
    d <- overlapped_int[overlapped_int$chr == cur_chr, ]
    cov <- calculate_coverage_for_overlapped(
      data=d, int_start=unique_combs[i, ]$interval_start, 
      int_end=unique_combs[i, ]$interval_end, win_len = win_len)
    res <- data.frame("chr" = cur_chr, "from" = unique_combs[i, ]$interval_start,
                      "to"=unique_combs[i, ]$interval_end, "cov"=cov, 
                      stringsAsFactors = FALSE)
    return(res)
  }
  
  data_overlapped <- data_overlapped %>%
    select(chr, from, to, cov) %>% 
    setNames(c("chr", "from", "to", feat_name))
  
  all_data <- rbind(data_grouped, data_overlapped) 
  
  return(all_data)
}



# win_len <- 10000
win_len <- 500000

output_path <- "data/features/"

path_to_data_folder <- "../cbp_data/data/"

mapping_files_nm <- list()
mapping_files_nm  <- c(
  'genome_regions/all_genome_regions.csv', 
  'DNA_methylation/preprocessed/data_all.csv',
  'histones/preprocessed/data_all.csv',
  "regulatory_regions/preprocessed/data_all.csv",
  "TAD_chromatin/all_tad_boundaries.csv",
  "transcr_factors/preprocessed/data_all.csv",
  "a-phased_repeats/data_all.csv",
  "direct_repeats/data_all.csv",
  "inverted_repeats/data_all.csv",
  "mirror_repeats/data_all.csv",
  "short_tandem_repeats/data_all.csv",
  "z-dna_motifs/data_all.csv",
  'G-quadruplexes/data_all.csv',
  'stemloops/sl_long_data_all.csv',
  'stemloops/sl_short_chr_1.csv',
  'stemloops/sl_short_chr_2.csv',
  'stemloops/sl_short_chr_3.csv',
  'stemloops/sl_short_chr_4.csv',
  'stemloops/sl_short_chr_5.csv',
  'stemloops/sl_short_chr_6.csv',
  'stemloops/sl_short_chr_7.csv',
  'stemloops/sl_short_chr_8.csv',
  'stemloops/sl_short_chr_9.csv',
  'stemloops/sl_short_chr_10.csv',
  'stemloops/sl_short_chr_11.csv',
  'stemloops/sl_short_chr_12.csv',
  'stemloops/sl_short_chr_13.csv',
  'stemloops/sl_short_chr_14.csv',
  'stemloops/sl_short_chr_15.csv',
  'stemloops/sl_short_chr_16.csv',
  'stemloops/sl_short_chr_17.csv',
  'stemloops/sl_short_chr_18.csv',
  'stemloops/sl_short_chr_19.csv',
  'stemloops/sl_short_chr_20.csv',
  'stemloops/sl_short_chr_21.csv',
  'stemloops/sl_short_chr_22.csv',
  'stemloops/sl_short_chr_X.csv'
  )

names(mapping_files_nm) <- c(
  "genome_regions",
  "DNA_methylation",
  "histones",
  "chromatin_state",
  "tad",
  "TF",
  "a-phased_repeats",
  "direct_repeats",
  "inverted_repeats",
  "mirror_repeats",
  "short_tandem_repeats",
  "z-dna_motifs",
  "G_quadruplexes",
  "sl_long",
  "sl_short_chr_1",
  "sl_short_chr_2",
  "sl_short_chr_3",
  "sl_short_chr_4",
  "sl_short_chr_5",
  "sl_short_chr_6",
  "sl_short_chr_7",
  "sl_short_chr_8",
  "sl_short_chr_9",
  "sl_short_chr_10",
  "sl_short_chr_11",
  "sl_short_chr_12",
  "sl_short_chr_13",
  "sl_short_chr_14",
  "sl_short_chr_15",
  "sl_short_chr_16",
  "sl_short_chr_17",
  "sl_short_chr_18",
  "sl_short_chr_19",
  "sl_short_chr_20",
  "sl_short_chr_21",
  "sl_short_chr_22",
  "sl_short_chr_X"
  )

chrs <-
  c(
    "1","2","3","4","5","6","7","8","9","10",
    "11","12","13","14","15","16","17","18","19","20","21","22","X"
  )


bins_path <- "data/good_bins/"
bins <- read.csv(paste0(bins_path, "bins_", format(win_len, scientific = FALSE), ".csv"), stringsAsFactors = FALSE, 
                 row.names = 1)


#######################################
# There are two different pipelines:
## 1. For conserved features (features which are same for all tissues biologically or technically)
## 2. For tissue-specific features

conserved_features <- c(
                        # "genome_regions",
                        # "tad",
                        # "G_quadruplexes",
                        # "sl_long",
                        # "z-dna_motifs",
                        # "a-phased_repeats",
                        # "direct_repeats",
                        # "inverted_repeats",
                        "mirror_repeats",
                        "short_tandem_repeats",
                        "sl_short_chr_1",
                        "sl_short_chr_2",
                        "sl_short_chr_3",
                        "sl_short_chr_4",
                        "sl_short_chr_5",
                        "sl_short_chr_6",
                        "sl_short_chr_7",
                        "sl_short_chr_8",
                        "sl_short_chr_9",
                        "sl_short_chr_10",
                        "sl_short_chr_11",
                        "sl_short_chr_12",
                        "sl_short_chr_13",
                        "sl_short_chr_14",
                        "sl_short_chr_15",
                        "sl_short_chr_16",
                        "sl_short_chr_17",
                        "sl_short_chr_18",
                        "sl_short_chr_19",
                        "sl_short_chr_20",
                        "sl_short_chr_21",
                        "sl_short_chr_22",
                        "sl_short_chr_X"
                        )

tissue_specific_features <- c("DNA_methylation",
                              "histones",
                              "chromatin_state",
                              "TF")


######################### CONSERVED FEATURES

total <- length(conserved_features)
pb <- tkProgressBar(
  title = "progress bar",
  min = 0,
  max = total,
  width = 300
)
i <- 1

for (feat in conserved_features){
  
  print(feat)
  
  # get data
  fl <- mapping_files_nm[feat]
  full_path <- paste0(path_to_data_folder, fl)
  data <- read.csv(file=full_path, header = TRUE, stringsAsFactors = FALSE)
  if ("X" %in% names(data)){
    data$X <- NULL
  }
  if ("V1" %in% names(data)){
    data$V1 <- NULL
  }
  data$chr <- as.character(data$chr)
  data <- data[data$chr %in% chrs, ]
  
  # create folder for specific features if it doesnot exist
  new_output_path <- paste0(output_path, feat)
  ifelse(!dir.exists(new_output_path), dir.create(new_output_path), FALSE)
  new_output_path <- paste0(new_output_path, "/")
  
  if (!'feature' %in% names(data)){
    data$feature <- feat
  }
  
  for (feat_name in unique(data$feature)){
    data_feature <- data[data$feature == feat_name, ]
    
    # get intersection with good bins
    data_feature$feature <- NULL
    data_intersected <- get_intersection_intervals(bins, data_feature)
    names(data_intersected) <- c('chr', 'from', 'to', 'chr_1', 'start', 'end')
    
    data_intersected <- data_intersected %>%
      mutate(
        corr_start = purrr::pmap_dbl(list(start, from), max),
        corr_end = purrr::pmap_dbl(list(end, to), min)
      ) %>% 
      arrange(chr, from, start)
    
    # check if features are overlapping
    
    res <- data_intersected %>%
      group_by(chr) %>%
      mutate(
        lag_end = lag(corr_end),
        lag_chr = lag(chr),
        lag_from = lag(from)
      ) %>%
      mutate(diff = corr_start - lag_end) %>%
      filter((chr == lag_chr) & (diff < 0) & (from == lag_from))
    
    if (nrow(res) == 0){
      print('no overlappping intervals')
      all_data <- get_coverage_not_overlapped(data_intersected, win_len, feat_name)
    } else {
      print('overlappping intervals')
      all_data <- get_coverage_overlapped(res, data_intersected, win_len, feat_name)    
    }
    
    write.csv(all_data, file = paste0(new_output_path, feat_name, "_", 
                                      format(win_len, scientific = FALSE), ".csv"))
    
  }
  
  setTkProgressBar(pb, i, label = paste(round(i / total * 100, 0), "% done"))
  i <- i + 1
}
close(pb)





######################### TSSUE-SPECIFIC FEATURES


feat <- tissue_specific_features[2]

# read data
fl <- mapping_files_nm[feat]
full_path <- paste0(path_to_data_folder, fl)
data <- read.csv(file=full_path, header = TRUE, stringsAsFactors = FALSE)
if ("X" %in% names(data)){
  data$X <- NULL
}
if ("V1" %in% names(data)){
  data$V1 <- NULL
}
data$chr <- as.character(data$chr)
data <- data[data$chr %in% chrs, ]

# create folder for specific features
new_output_path <- paste0(output_path, feat)
ifelse(!dir.exists(new_output_path), dir.create(new_output_path), FALSE)
new_output_path <- paste0(new_output_path, "/")

if (!'feature' %in% names(data)){
  data$feature <- feat
}

for (feat_name in unique(data$feature)){

  data_feature <- data[data$feature == feat_name, ]
  
  # get intersection with good bins
  data_feature$feature <- NULL
  
  cancer_types <- unique(data_feature$cancer_type)
  all_cancer_data <- data.frame()
  
  for (cur_cancer_type in cancer_types){
    
    data_feature_cancer <- data_feature[data_feature$cancer_type == cur_cancer_type, ]
    data_feature_cancer$cancer_type <- NULL
    # get intersection with good bins
    data_feature_cancer$feature <- NULL
    data_intersected <- get_intersection_intervals(bins, data_feature_cancer)
    names(data_intersected) <- c('chr', 'from', 'to', 'chr_1', 'start', 'end')
    
    data_intersected <- data_intersected %>%
      mutate(
        corr_start = purrr::pmap_dbl(list(start, from), max),
        corr_end = purrr::pmap_dbl(list(end, to), min)
      ) %>% 
      arrange(chr, from, start)
    
    # check if features are overlapping
    
    res <- data_intersected %>%
      group_by(chr) %>%
      mutate(
        lag_end = lag(corr_end),
        lag_chr = lag(chr),
        lag_from = lag(from)
      ) %>%
      mutate(diff = corr_start - lag_end) %>%
      filter((chr == lag_chr) & (diff < 0) & (from == lag_from))
    
    if (nrow(res) == 0){
      print('no overlappping intervals')
      all_data <- get_coverage_not_overlapped(data_intersected, win_len, feat_name)
    } else {
      print('overlappping intervals')
      all_data <- get_coverage_overlapped(res, data_intersected, win_len, feat_name)    
    }
    
    all_data$cancer_type <- cur_cancer_type
    
    all_cancer_data <- rbind(all_cancer_data, all_data)
  }
  
  write.csv(all_cancer_data, file = paste0(new_output_path, feat_name, "_", win_len, ".csv"))

}
