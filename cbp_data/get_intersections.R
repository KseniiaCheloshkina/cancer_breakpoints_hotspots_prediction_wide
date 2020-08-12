library(dplyr)
library(reshape2)
require(data.table)

# function to get the number of points from df_points 
# located in intervals from df_intervals
# df_points: data.table(chr, chr_bkpt, cancer_type, id)
# df_intervals: data.table(chr, start, end) 

get_intersection_one_point <- function(df_points, df_intervals){
  
  joined_df <- df_points[df_intervals, 
                         on=.( 
                           chr_bkpt <= end, 
                           chr_bkpt >= start),
                         allow.cartesian=TRUE
                         ] %>%
    filter(!is.na(id)) %>%
    select(chr, cancer_type, id) %>%
    group_by(chr, cancer_type) %>%
    summarize(n_bkpt = n_distinct(id))
  
  return(data.frame(joined_df))
}




# function to get intervals intersections of two dataframes   
# INPUT:
# df_intervals_1: data.frame(chr, start, end)
# df_intervals_2: data.frame(chr, start, end) 
# RETURN:
# df_result: data.frame(chr, start, end, chr_1,start_1, end_1)
# where each row is intersected intervals of df_intervals_1 and df_intervals_2

get_intersection_intervals <- function(df_intervals_1, df_intervals_2){
  
  rownames(df_intervals_1) <- NULL
  rownames(df_intervals_2) <- NULL
  df_intervals_1$id <- rownames(df_intervals_1)
  df_intervals_2$id <- rownames(df_intervals_2)
  df_intervals_2 <- df_intervals_2 %>% 
    rename(chr_1 = chr, start_1 = start, end_1 = end, id_1 = id)
  
  setDT(df_intervals_1)
  setDT(df_intervals_2)
  chr_list <- unique(df_intervals_1[, chr])
  
  merged_data <- data.frame()
  
  for (chr in chr_list){
    df_part_1 <- df_intervals_1[chr, on="chr"]
    df_part_2 <- df_intervals_2[chr, on="chr_1"]
    
    joined_df_1_s <- df_part_1[df_part_2, 
                               on=.( start <= start_1, 
                                     end > start_1),
                         allow.cartesian=TRUE] %>%
      filter(!is.na(id)) %>%
      filter(!is.na(id_1))  %>%
      select(id, id_1)
    
    joined_df_1_e <- df_part_1[df_part_2, 
                               on=.( start < end_1, 
                                     end >= end_1),
                         allow.cartesian=TRUE]  %>%
      filter(!is.na(id)) %>%
      filter(!is.na(id_1))  %>%
      select(id, id_1)
    
    joined_df_2_s <- df_part_1[df_part_2, 
                               on=.( start >= start_1, 
                                     start < end_1),
                         allow.cartesian=TRUE] %>%
      filter(!is.na(id)) %>%
      filter(!is.na(id_1))  %>%
      select(id, id_1)
    
    final_df <- rbind(joined_df_1_s, joined_df_1_e, joined_df_2_s) 
    final_df <- data.frame(final_df)
    final_df <- unique(final_df)
    
    final_df1 <- final_df %>%
      inner_join(df_part_1, by = "id") %>%
      inner_join(df_part_2, by = "id_1") %>%
      select(-c(id, id_1))
    
    merged_data <- rbind(merged_data, final_df1)
    
  }
  
  return(merged_data)
  
}
