script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(ggplot2)
library(reshape2)

source("../cbp_data/get_intersections.R")


# set window size
win_len <- 10000

path_to_data_folder <- "../cbp_data/data/"
path_to_bkpt <- 'raw breakpoints/raw_release_28/after_preprocessing/'

output_path <- "data/target/"

# read bins
bins_path <- "data/good_bins/"
bins <- read.csv(paste0(bins_path, "bins_", format(win_len, scientific = FALSE), ".csv"), 
                 stringsAsFactors = FALSE, row.names = 1)

# read breakpoints data
cancer_types <- list.files(paste0(path_to_data_folder, path_to_bkpt))
cancer_types <- cancer_types[grep("_all_data.csv", cancer_types)]

bkpt_data <- data.frame()
for (i in cancer_types) {
  canc <- strsplit(i, "_")[[1]][1]
  data <- read.csv(paste0(path_to_data_folder, path_to_bkpt, i),
                   stringsAsFactors = FALSE)
  data <- data %>%
    select(chr, chr_bkpt_beg, chr_bkpt_end)
  
  data$cancer_type <- canc
  bkpt_data <- rbind(bkpt_data, data)
}

bkpt_data <- bkpt_data %>%
  setNames(c("chr", "start", "end", "cancer_type"))

# get intersection with good bins
# for each cancer type
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
unique(all_bkpt_counts$chr)
all_bkpt_counts <- all_bkpt_counts[all_bkpt_counts$chr != "Y", ]

write.csv(all_bkpt_counts, file = paste0(output_path, "all_data_", win_len, ".csv"), 
          row.names = FALSE)
          