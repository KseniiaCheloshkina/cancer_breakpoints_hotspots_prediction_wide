Genome-wide study of cancer breakpoints hotspots

Structure of repository:
- run: all scripts used in paper. Main scripts:
  - filter_bins.R: script for getting the list of all genome windows for specified window length excluding first and last windows, centromeres, telomeres and blacklisted regions
  - target_<N_name>.R: scripts for processing breakpoints data from the ICGC (filtering, merging, getting breakpoints density and hotspots)
  - calculate_coverage.R: script for calculating coverage of features in genome windows (has two pipelines: for "conserved" and tissue-specific features)
  - create_dataset.R: script for collecting targets and features data into one single file by window length
  - EDA dataset.Rmd: RMarkdown with exploratory data analysis of dataset
  - train.R: script for binary classification model selection and checking features contribution

This repository has [KseniiaCheloshkina/cbp_data](https://github.com/KseniiaCheloshkina/cbp_data/) repository as a submodule, using some general purpose tools from its run/ directory as well as all raw data downloaded from sources from data/ folder.
