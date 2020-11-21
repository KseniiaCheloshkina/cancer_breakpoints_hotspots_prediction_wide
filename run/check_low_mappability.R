library(dplyr)
library(reshape2)
library(stringr)
require(data.table)

script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)

source("../cbp_data/get_intersections.R")

win_len <- 100000

## load bins (excluded dac blacklist and first/last bin of each chromosomes,
# telomeres and centromeres)
bins <- read.csv(file="../data/good_bins/bins_100000.csv")
bins$X <- NULL
bins$chr <- as.character(bins$chr)
bins <- bins[bins$chr != "Y", ]

## load target
df_target <- read.csv(file="../data/target/all_data_10000.csv")
df_target <- df_target[c("chr", "from", "to", "bkpt_in_window_breast")]

bkpt_data <- df_target
bkpt_data$new_from <- floor(bkpt_data$from / win_len) * win_len
bkpt_data$new_to <- bkpt_data$new_from + win_len
bkpt_cols <- grep(x = names(bkpt_data), pattern = "bkpt", value = TRUE)
bkpt_data <- bkpt_data %>%
  group_by(chr, new_from, new_to) %>%
  summarise_at(bkpt_cols, sum)
bkpt_data <- bkpt_data %>%
  rename("end" = "new_to") %>%
  rename("start" = "new_from")

bins <- bins %>%
  left_join(
    bkpt_data, by=c("chr", "start", "end")
  )


## load blacklist
low_map_path <-'../../cbp_data/data/genome/low_mappability.bedGraph'
df_low_map <- read.table(low_map_path, header=FALSE)
names(df_low_map) <- c("chr", "start", "end", "score")

# filter Y chromosome
df_low_map$chr <- as.character(df_low_map$chr)
df_low_map <- df_low_map[!df_low_map$chr %in% c("chrM", "chrY"), ]
df_low_map$chr <- sapply(df_low_map$chr, function(x) gsub(pattern = "chr", replacement = "", x=x))


# The Duke Uniqueness track 35mer
# Scores are normalized to between 0 and 1, with 1 representing a 
# completely unique sequence and 0 representing a sequence that occurs
# more than 4 times in the genome (excluding chrN_random and alternative 
# haplotypes). A score of 0.5 indicates the sequence occurs exactly twice,
# likewise 0.33 for three times and 0.25 for four times.

# select regions with uniqueness score lower than threshold 'thr'
# these regions are called 'low-mappable' (are met more frequently than threshold)
# and should be removed

thrs <- c(0.25, 0.33, 0.5)
all_data <- data.frame()

for (thr in thrs){
  df_low_map_thr <- df_low_map[df_low_map$score < thr, ]
  df_low_map_thr$score <- NULL
  
  windows_with_blacklisted <- get_intersection_intervals(bins, df_low_map_thr)
  windows_with_blacklisted <- windows_with_blacklisted %>%
    mutate(
      start_final = ifelse(start_1 < start, start, start_1),
      end_final = ifelse(end_1 > end, end, end_1)
    )
  windows_with_blacklisted <- windows_with_blacklisted %>%
    mutate(len = end_final - start_final) %>%
    group_by(chr, start, end) %>%
    summarize(
      len_low_mappable = sum(len)
    ) %>%
    mutate(
      len_total = end-start,
      coverage_low_mappable = len_low_mappable / len_total
    ) %>%
    select(-len_total)
  
  df_final <- bins %>%
    left_join(windows_with_blacklisted) %>%
    mutate(
      thr = thr
    )
  all_data <- rbind.data.frame(all_data, df_final)
  write.csv(all_data, "../data/output_third/duke_low_mappable_all.csv")
}



# by chr
chrs <- unique(df_low_map$chr)
thr <- 1
df_low_map <- df_low_map[df_low_map$score < thr, ]
all_data <- data.frame()
for (chr in chrs){
  print(chr)
  df_low_map_thr <- df_low_map[df_low_map$chr == chr, ]
  df_low_map_thr$score <- NULL
  bins_cur <- bins[bins$chr == chr, ]
  windows_with_blacklisted <- get_intersection_intervals(bins_cur, df_low_map_thr)
  windows_with_blacklisted <- windows_with_blacklisted %>%
    mutate(
      start_final = ifelse(start_1 < start, start, start_1),
      end_final = ifelse(end_1 > end, end, end_1)
    )
  windows_with_blacklisted <- windows_with_blacklisted %>%
    mutate(len = end_final - start_final) %>%
    group_by(chr, start, end) %>%
    summarize(
      len_low_mappable = sum(len)
    ) %>%
    mutate(
      len_total = end-start,
      coverage_low_mappable = len_low_mappable / len_total
    ) %>%
    select(-len_total)
  
  df_final <- bins_cur %>%
    left_join(windows_with_blacklisted) %>%
    mutate(
      thr = thr
    )
  all_data <- rbind.data.frame(all_data, df_final)
  write.csv(all_data, "../data/output_third/duke_low_mappable_all_1.csv")
}



## analytics
library(ggplot2)
all_data <- read.csv("../data/output_third/duke_low_mappable_all.csv")
all_data$X <- NULL

score_thrs <- c(0.25, 0.33, 0.5)
bins_ratio_unmappable_thr <- c(0.01, 0.05, 0.1, 0.15, 0.2)
all_res <- data.frame()
for (thr_score in score_thrs){
  data_part <- all_data %>%
    filter(
      thr == thr_score
    )
  for (thr_bins in bins_ratio_unmappable_thr){
    thr_bins_summary <- data_part %>%
      filter(
        coverage_low_mappable < thr_bins
      ) %>%
      group_by(thr) %>%
      summarize(
        n_bins = n(),
        n_breaks = sum(bkpt_in_window_breast)
      ) %>%
      mutate(
        thr_cov_bins = thr_bins
      )
    all_res <- rbind.data.frame(all_res, thr_bins_summary)
  }
}

data_part_summary <- data_part %>%
  summarize(
    total_n_bins = n(),
    total_n_breaks = sum(bkpt_in_window_breast)
  )

all_res <- all_res %>%
  mutate(
    total_n_bins = data_part_summary$total_n_bins,
    total_n_breaks = data_part_summary$total_n_breaks,
    proportion_breaks_retained = n_breaks / total_n_breaks,
    proportion_bins_retained = n_bins / total_n_bins
  )
all_res$thr_cov_bins <- as.character(all_res$thr_cov_bins)
all_res$proportion_breaks_retained <- round(all_res$proportion_breaks_retained, 3)
all_res$proportion_bins_retained <- round(all_res$proportion_bins_retained, 3)
map_alias_thr <- data.frame(thr = c(0.25, 0.33, 0.5), thr_names = c("more than 4 times",
                                                                     "more than 3 times",
                                                                     "more than 2 times"))
all_res <- all_res %>% inner_join(map_alias_thr)

ggplot(all_res, aes(x=thr_names, y=thr_cov_bins, fill=proportion_breaks_retained)) +
  geom_tile()+
  geom_text(aes(label=proportion_breaks_retained))+
  xlab("Exclude regions with uniqueness")+
  ylab("Coverage of excluded regions less than")

ggplot(all_res, aes(x=thr_names, y=thr_cov_bins, fill=proportion_bins_retained)) +
  geom_tile()+
  geom_text(aes(label=proportion_bins_retained))+
  xlab("Exclude regions with uniqueness")+
  ylab("Coverage of excluded regions less than")

