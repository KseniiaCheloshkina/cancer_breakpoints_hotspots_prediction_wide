This project contains source code for the paper of Cheloshkina K., Poptsova M. "Comprehensive analysis of cancer breakpoint hotspots with machine learning models". The main goal of the study is to estimate prediction power of different genomic and epigenomic features for prediction of cancer breakpoints hotspots.

Structure of repository:
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
