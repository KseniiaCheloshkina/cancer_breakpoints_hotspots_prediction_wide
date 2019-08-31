library(dplyr)
library(reshape2)
library(stringr)
require(data.table)

script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)

source("../cbp_data/get_intersections.R")
source("../cbp_data/get_genome_windows.R")

# set window size
win_len <- 1000000

path_to_data_folder <- "../../cbp_data/data/"

centromeres_path <- paste0(path_to_data_folder, "genome/gaps.bed")
blacklist_path <- paste0(path_to_data_folder, "genome/dac_blacklist.bed")

output_path <- "../data/good_bins/"

############### Regions to filter

# get genome windows for specified window length 
chr_windows <- get_chr_windows_hg19(win_len)

# remove first and last bins of ech chromosome
chr_windows <- chr_windows %>%
  group_by(chr) %>%
  mutate(
    win_rank_asc = dense_rank(from),
    win_rank_desc = dense_rank(desc(from))
  ) %>%
  ungroup() %>%
  filter(
    (win_rank_asc != 1) &
    (win_rank_desc != 1)
  ) %>%
  select(chr, from, to)

names(chr_windows)<- c('chr', 'start', 'end')

# get intersections with centromeres and telomeres
gaps <- read.table(centromeres_path,stringsAsFactors = FALSE)
gaps <- gaps[gaps$V8 %in% c('centromere', 'telomere'), ]
gaps <- gaps[, 2:4]
names(gaps) <- c('chr', 'start', 'end')
gaps$chr <- gsub(x = gaps$chr, pattern = "chr", replacement = "")

windows_with_gaps <- get_intersection_intervals(chr_windows, gaps)
windows_with_gaps$is_intersected <- 1
windows_with_gaps <- windows_with_gaps %>%
  select(chr, start, end, is_intersected)

chr_windows <- chr_windows %>% 
  left_join(windows_with_gaps, by=c("chr", 'start', 'end')) %>% 
  filter(is.na(is_intersected)) %>% 
  select(-is_intersected)


# get intersections with blacklisted regions
dac_bl <- read.table(blacklist_path, header=FALSE,stringsAsFactors = FALSE)
dac_bl <- dac_bl[2:4]
names(dac_bl) <- c('chr', 'start', 'end')
dac_bl$chr <- sapply(dac_bl$chr, function(x) gsub(pattern = "chr", replacement = "", x=x))

windows_with_blacklisted <- get_intersection_intervals(chr_windows, dac_bl)
windows_with_blacklisted <- windows_with_blacklisted %>%
  mutate(is_intersected = 1)  %>%
  select(chr, start, end, is_intersected)

chr_windows <- chr_windows %>% 
  left_join(windows_with_blacklisted, by=c("chr", 'start', 'end')) %>% 
  filter(is.na(is_intersected)) %>% 
  select(-is_intersected)

write.csv(chr_windows, paste0(output_path, "bins_", format(win_len, scientific = FALSE),".csv"))
