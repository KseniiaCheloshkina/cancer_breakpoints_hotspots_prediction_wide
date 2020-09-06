script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd('..')

library(ggplot2)
library(dplyr)
source("../cbp_data/get_intersections.R")

#### READ GENES DATA

df <- read.table('data/adhoc/genes.txt')
names(df) <- c("chr", "strand", "start", "end", "name")

# select for each gene max end and then max start
df <- df %>%
  group_by(name) %>%
  mutate(
    end_rank = dense_rank(desc(end))
  ) %>%
  filter(
    end_rank == 1
  ) %>%
  unique() %>%  
  mutate(
    start_rank = dense_rank(desc(start)),
  ) %>%
  filter(
    start_rank == 1
  ) %>%
  group_by(name) %>%
  mutate(rn = row_number()) %>%
  filter(rn == 1)

# correct with upstream
df[df$strand == "+", "start"] <- df[df$strand == "+", "start"] - 1000
df[df$strand == "-", "end"] <- df[df$strand == "-", "end"] + 1000
df[c("strand", "end_rank", "start_rank", "rn")] <- NULL

# select chromosomes
df$chr <- gsub(x = df$chr, pattern = "chr", replacement = "")
good_chr <- c(as.character(seq(1, 22)), "X")
df <- df[df$chr %in% good_chr, ]


#### READ BREAKPOINTS DATA

bkpt <- read.csv("data/datasets/dataset_100000.csv")
ss_cols <- c("chr", "from", "to")
hsp_cols <- grep(names(bkpt), pattern = "hsp", value=T)
needed_cols <- c(ss_cols, hsp_cols)
bkpt <- bkpt[needed_cols]
bkpt <- bkpt %>%
  rename(start = from, end = to)
ss_cols <- c("chr", "start", "end")
bkpt$chr <- as.character(bkpt$chr)



#### intersect with breakpoints and hotspots

all_df <- data.frame()

for (target_col in hsp_cols){
  cancer_type <- strsplit(target_col, "_")[[1]][3]
  agg_level <- strsplit(target_col, "_")[[1]][2]
  bkpt_part <- bkpt[bkpt[target_col] == 1, ss_cols]
  bkpt_genes <- get_intersection_intervals(bkpt_part, df)
  genes <- bkpt_genes %>%
    select(name) %>%
    unique() %>%
    mutate(
      cancer_type = cancer_type,
      agg_level = agg_level,
      gene_name = as.character(name)
      ) %>%
    select(-name)
  all_df <- rbind.data.frame(all_df, genes)
}
write.csv(all_df, "data/adhoc/genes_intersection.csv")


all_df %>%
  group_by(gene_name) %>%
  summarize(n = n()) %>%
  arrange(desc(n)) %>%
  first(100) %>%
  ggplot(aes(x=gene_name, y=n)) + geom_bar(stat='identity')


all_df <- read.csv("data/adhoc/genes_intersection.csv")
unique_df <- all_df %>% 
  select(agg_level, gene_name) %>%
  unique()

unique_df %>%
  group_by(agg_level) %>%
  summarize(n())

write.csv(unique_df, "data/adhoc/genes_intersection_unique_all_cancers.csv")
