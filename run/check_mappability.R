library(dplyr)
library(reshape2)
library(stringr)
require(data.table)

script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)

#### download file 
# 
#### download utility to convert bigWig binary file to bed
# rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bigWigToBedGraph ./
# cd /
# chmod +x /home/islam/bigWigToBedGraph
##### convert
# /home/islam/bigWigToBedGraph /home/islam/cbp_data/data/genome/wgEncodeDukeMapabilityUniqueness35bp.bigWig /home/islam/out.bedGraph



# get intersections with blacklisted regions
dac_bl <- read.table("../../cbp_data/data/genome/wgEncodeDukeMapabilityUniqueness35bp.bigWig", nrows = 2,
                     stringsAsFactors = FALSE)
dac_bl %>% head()
dac_bl <- dac_bl[2:4]
names(dac_bl) <- c('chr', 'start', 'end')
dac_bl$chr <- sapply(dac_bl$chr, function(x) gsub(pattern = "chr", replacement = "", x=x))
