script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(dplyr)
library(openxlsx)

output_folder <- "data/target/"

# set window size
win_len <- 100000

# read breakpoints data
path_to_data_folder <- "data/target/all_data_10000.csv"
data <- read.csv(path_to_data_folder, stringsAsFactors = F)
select_cols <- c("chr", "from", "to",  
                 "bkpt_in_window_blood", "bkpt_in_window_bone", "bkpt_in_window_brain",
                 "bkpt_in_window_breast", "bkpt_in_window_liver", "bkpt_in_window_ovary",
                 "bkpt_in_window_pancreatic", "bkpt_in_window_prostate", 
                 "bkpt_in_window_skin", "bkpt_in_window_uterus")
data <- data[select_cols]

# aggregate to win_len
if (win_len == 10000){
  bkpt_data <- data
} else {
  bkpt_data <- data
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



# discover quantiles of distribution
bkpt_cols <- grep(x = select_cols, "bkpt", value=T)
q_borders <- data.frame()

for (bkpt_col in bkpt_cols){
  
  cancer_type <- strsplit(bkpt_col, "_")[[1]][[4]]
  data_part <- bkpt_data[bkpt_data[bkpt_col] > 0, ]
  s <- data.frame(t(data.frame(quantile(data_part[, bkpt_col], c(0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 0.95)))))
  rownames(s) <- NULL
  s$cancer_type <- cancer_type
  q_borders <- rbind(q_borders, s)
}


# select specific quantiles
q_for_cancer <- list()

q_part <- list(2, 3, 4)
names(q_part) <- c("75.", "90.", "95.")
q_for_cancer[['blood']] <- q_part

q_part <- list(2, 3, 4)
names(q_part) <- c("50.", "90.", "95.")
q_for_cancer[['bone']] <- q_part

q_part <- list(2, 4)
names(q_part) <- c("75.", "95.")
q_for_cancer[['brain']] <- q_part

q_part <- list(5, 9, 14, 19)
names(q_part) <- c("50.", "75.", "90.", "95.")
q_for_cancer[['breast']] <- q_part

q_part <- list(2, 4, 5)
names(q_part) <- c("50.", "90.", "95.")
q_for_cancer[['liver']] <- q_part

q_part <- list(3, 5, 7, 9)
names(q_part) <- c("50.", "75.", "90.", "95.")
q_for_cancer[['ovary']] <- q_part

q_part <- list(3, 5, 8, 11)
names(q_part) <- c("50.", "75.", "90.", "95.")
q_for_cancer[['pancreatic']] <- q_part

q_part <- list(5, 9, 13,17)
names(q_part) <- c("50.", "75.", "90.", "95.")
q_for_cancer[['prostate']] <- q_part

q_part <- list(2,4,6,9)
names(q_part) <- c("50.", "75.", "90.", "95.")
q_for_cancer[['skin']] <- q_part

q_part <- list(2, 4, 6)
names(q_part) <- c("50.", "90.", "95.")
q_for_cancer[['uterus']] <- q_part


# get labels for them and create cols
for (bkpt_col in bkpt_cols){
  
  cancer_type <- strsplit(bkpt_col, "_")[[1]][[4]]
  q_data <- q_for_cancer[[cancer_type]]
  for (q in names(q_data)){
    new_col <- paste0("q_", q, "_", cancer_type)
    bkpt_data[new_col] <- ifelse(bkpt_data[bkpt_col] >= q_data[q], 1, 0)
  }
}

bkpt_data_new <- bkpt_data %>% select(-bkpt_cols)

# save stats
stats <- bkpt_data_new %>%
  select(-c("chr","from", "to")) %>%
  summarize_all(sum) %>%
  t %>%
  as.data.frame()
names(stats) <- "number of positive examples"
stats$cancer_type <- unlist(lapply(strsplit(row.names(stats), "_"), function(x) x[[2]]))
stats$labeling <- unlist(lapply(strsplit(row.names(stats), "_"), function(x) x[[3]]))
wb <- createWorkbook()
addWorksheet(wb, "num_win")
writeData(wb, sheet="num_win", stats)
saveWorkbook(wb, "reports/bkpt_randomness.xlsx", overwrite = T)

# save data
write.csv(bkpt_data_new, file=paste0(output_folder, "q_bkpt_", format(win_len, scientific = F),
                                 ".csv"), row.names = F)


#### get target "hotspots vs breakpoints"

# read final target data
path_to_data_folder <- "data/datasets/dataset_100000.csv"
data <- read.csv(path_to_data_folder, stringsAsFactors = F)
hsp_cols <- grep(names(data), pattern = "hsp", value = T)
data <- data %>%
  select(c("chr", "from", "to", hsp_cols))

# merge with number of breakpoints in window
bkpt_data_bkpt <- bkpt_data %>% 
  select(c("chr", "from", "to", bkpt_cols))
all_data <- data %>%
  inner_join(
    bkpt_data_bkpt, by=c("chr","from", "to")
  )
# only for hotspots 
hsp_cols <- hsp_cols[-grep(x = hsp_cols, pattern = 'all')]
# for each target create a new target col with labeled windows without breakpoints
for (col in hsp_cols){
  
  cancer_type <- strsplit(col, "_")[[1]][3]
  new_col <- paste0(col, "_red")
  all_data[new_col] <- all_data[col]
  
  bkpt_col <- grep(x = bkpt_cols, pattern = cancer_type, value=T)
  all_data[all_data[bkpt_col] == 0, new_col] <- 2
}


target_cols <- c("chr", "from", "to", grep(x = names(all_data), pattern="red", value=T)) 
write.csv(all_data[target_cols], file=paste0(output_folder, "hsp_bkpt_", format(win_len, scientific = F),
                                     ".csv"), row.names = F)

