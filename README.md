This project contains source code for the paper of Cheloshkina K., Poptsova M. "Comprehensive analysis of cancer breakpoints reveals signatures of genetic and epigenetic contribution to cancer genome rearrangements". The main goal of the study is to estimate prediction power of different genomic and epigenomic features for prediction of cancer breakpoints hotspots.
### Datasets
Data used for the study are publicly available. For the purpose of models benchmarking and making the research results reproducible we upload in the repository  genomic coordinates of train-test splits. In the study we used 30-times repeated train-test split withhelding 30% of sample in test set for each repeat.
The data are located in ` datasets/ ` folder:
 - `classification/` folder contains data for the task of classification of genome windows as breakpoints hotspot or not. 
 - `randomness/` folder contains data for the task of breakpoints randomness nature investigation. 
 The folder contains files of 2 types: "hsp_99.5._blood_red.csv" (files for classification of windows with breakpoints as hotspots or individual breakpoints) and "q_50._bone.csv" (files for classification of genome windows as breakpoints hotspot or not for lower labeling types). 

File names comprised of 2 parts: "hsp_**99.9.**_**liver**.csv", where:
- **99.9.** denotes labeling type ("99.9" corresponds to labeling where genome windows with density higher than 99.9% of breakpoints density distribution are labeled as hotspots; "all" labeling type corresponds to breakpoints labels);
- **liver** corresponds to cancer type.

Each file contains next fields: "chr" (chromosome name), "from" (genome window start coordinate in specified chromosome),  "to" (genome window end coordinate in specified chromosome),	"data_type" (used in train or test),	"target" (target label), 	"cancer_type" (type of cancer),	"labeling_type" (ratio of all windows denoted as hotspots),	"repeat" (index of split in 30-times train-test resampling scheme).

  
### Structure of repository:
- run: all scripts used in paper. Main scripts:
  - filter_bins.R: script for getting the list of all genome windows for specified window length excluding first and last windows, centromeres, telomeres and blacklisted regions
  - target_<N_name>.R: scripts for processing breakpoints data from the ICGC (filtering, merging, getting breakpoints density and hotspots)
  - calculate_coverage.R: script for calculating coverage of features in genome windows (has two pipelines: for "conserved" and tissue-specific features)
  - create_dataset.R: script for collecting targets and features data into one single file by window length
  - model_selection.R: script for model selection (and hyperparameters optimization)
  - features.R: tools for feature preprocessing
  - check_features.R: script for feature engineering 
  - train.R: script for modeling using all features
  - boruta_fs.R: script for Boruta Feature Selection procedure execution 
  - train_less_features.R: script for modeling using specified feature set
  - check_significant_features.R: script for stepwise feature addition to Boruta feature set
  - train_by_feature_group.R: script for modeling using each feature group separately
  - systematic_target.R: script for creating new targets for checking the hypothesys of breakpoints randomness
  - train_bkpt_randomness.R: script for modeling using targets for checking the hypothesys of breakpoints randomness

This repository has [KseniiaCheloshkina/cbp_data](https://github.com/KseniiaCheloshkina/cbp_data/) repository as a submodule, using some general purpose tools from its run/ directory.
