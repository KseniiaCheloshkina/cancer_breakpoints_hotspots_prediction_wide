script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

data_folder <- "data/features/"


chr_list <- c(
  "DNA_methylation_chr_1",
  "DNA_methylation_chr_2",
  "DNA_methylation_chr_3",
  "DNA_methylation_chr_4",
  "DNA_methylation_chr_5",
  "DNA_methylation_chr_6",
  "DNA_methylation_chr_7",
  "DNA_methylation_chr_8",
  "DNA_methylation_chr_9",
  "DNA_methylation_chr_10",
  "DNA_methylation_chr_11",
  "DNA_methylation_chr_12",
  "DNA_methylation_chr_13",
  "DNA_methylation_chr_14",
  "DNA_methylation_chr_15",
  "DNA_methylation_chr_16",
  "DNA_methylation_chr_17",
  "DNA_methylation_chr_18",
  "DNA_methylation_chr_19",
  "DNA_methylation_chr_20",
  "DNA_methylation_chr_21",
  "DNA_methylation_chr_22",
  "DNA_methylation_chr_X"
)

chr_list <- c(
  "chromatin_state_chr_1",
  "chromatin_state_chr_2",
  "chromatin_state_chr_3",
  "chromatin_state_chr_4",
  "chromatin_state_chr_5",
  "chromatin_state_chr_6",
  "chromatin_state_chr_7",
  "chromatin_state_chr_8",
  "chromatin_state_chr_9",
  "chromatin_state_chr_10",
  "chromatin_state_chr_11",
  "chromatin_state_chr_12",
  "chromatin_state_chr_13",
  "chromatin_state_chr_14",
  "chromatin_state_chr_15",
  "chromatin_state_chr_16",
  "chromatin_state_chr_17",
  "chromatin_state_chr_18",
  "chromatin_state_chr_19",
  "chromatin_state_chr_20",
  "chromatin_state_chr_21",
  "chromatin_state_chr_22",
  "chromatin_state_chr_X"
)

chr_list <- c(
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


save_folder <- "DNA_methylation/"
base_name <- "DNA_methylation"

save_folder <- "chromatin_state/"
base_name <- "DNase_seq"

save_folder <- "sl_short/"
base_name <- "stemloops_6_15"

all_win_len <- c(10000, 20000, 50000, 100000, 500000, 1000000)

for (win_len in all_win_len){
  
  all_data <- data.frame()
  
  for (chr_data_path in chr_list){
    
    fl <- list.files(paste0(data_folder, chr_data_path))
    fl <- grep(x = fl, pattern = paste0("_", format(win_len, scientific = FALSE), ".csv"), value = TRUE)
    
    data <- read.csv(paste0(data_folder, chr_data_path, "/", fl), row.names = 1)
    all_data <- rbind(all_data, data) 
  }
  
  write.csv(all_data, 
            file = paste0(data_folder, save_folder, base_name, "_", 
                          format(win_len, scientific = FALSE), ".csv"),
            row.names = FALSE)
}
