#!/bin/bash


cores="$2"
input_bins_folder="$1"/binning/unique_bins
output_results="$1"/metrics/prokka

# loading conda environment
conda activate prokka-env

source MuDoGeR/installation/config_db

unset LD_PRELOAD

metawrap annotate_bins -t $cores -o $output_results -b $input_bins_folder
