script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
setwd("..")

# UNION BREAKPOINTS FROM SEVERAL FILES

# data path
input_folder <- "../cbp_data/data/breakpoints/raw_release_28/"
setwd(input_folder)

cancer_types <- list.dirs()
cancer_types <- cancer_types[cancer_types != "."]
cancer_types <- cancer_types[cancer_types != "./after_preprocessing"]


for (i in cancer_types){
  
  print(i)

  cancer_type_files <- list.files(i)

  data_frame_cols <- c("icgc_donor_id",
                     "project_code",
                     "icgc_sample_id",
                     "submitted_sample_id",
                     "submitted_matched_sample_id",
                     "variant_type",
                     "sv_id",
                     "assembly_version",
                     "chr_from",                        
                     "chr_from_bkpt",                   
                     "chr_from_strand"  ,               
                     "chr_from_range"  ,                
                     "chr_to",                          
                     "chr_to_bkpt"   ,                 
                     "chr_to_strand"  ,                 
                     "chr_to_range" ,
                     "quality_score",                   
                     "gene_affected_by_bkpt_from" ,     
                     "gene_affected_by_bkpt_to"    ,    
                     "transcript_affected_by_bkpt_from",
                     "transcript_affected_by_bkpt_to"  
  )
  all_data <- as.data.frame(x=t(rep(0, length(data_frame_cols))))
  names(all_data) <- data_frame_cols
  all_data <- all_data[-1, ]

  for (j in cancer_type_files){

    donor_file <- paste(i, j, sep="/")

    donor_data <- read.table(donor_file, sep="\t", fill = TRUE, header=TRUE)
    donor_data <- subset(donor_data, select = data_frame_cols)
    print(unique(donor_data$assembly_version))
    all_data <- rbind(all_data, donor_data)
  }

  unique_rows <- unique(all_data)
  
  write.csv(unique_rows, file=paste(gsub("./", "", i), "_all_data.csv", sep=""))
}


# DATA FOR ALL CANCERS IN ONE FILE

cancer_files <- list.files()
cancer_files <- cancer_files[grep("_all_data.csv", cancer_files)]


data_frame_cols <- c(
  "icgc_donor_id",
  "project_code",
  "icgc_sample_id",
  "submitted_sample_id",
  "submitted_matched_sample_id",
  "variant_type",
  "sv_id",
  "chr_from",
  "chr_from_bkpt",
  "chr_from_strand"  ,
  "chr_from_range"  ,
  "chr_to",
  "chr_to_bkpt"   ,
  "chr_to_strand"  ,
  "chr_to_range" ,
  "quality_score",
  "gene_affected_by_bkpt_from" ,
  "gene_affected_by_bkpt_to"    ,
  "transcript_affected_by_bkpt_from",
  "transcript_affected_by_bkpt_to"
)

all_data <- as.data.frame(x = t(rep(0, length(data_frame_cols))))
names(all_data) <- data_frame_cols
all_data <- all_data[-1, ]


for (i in cancer_files) {
  data <- read.csv(i)
  data <- subset(data, select = data_frame_cols)
  all_data <- rbind(all_data, data)
}

write.csv(all_data, file = "all_cancer_data.csv")
