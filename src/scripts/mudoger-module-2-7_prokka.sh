#!/bin/bash


cores="$2"
input_bins_folder="$1"/binning/unique_bins
output_results="$1"/metrics/prokka

# loading conda environment
conda activate metawrap-env

source Mudoger/installation/config 

unset LD_PRELOAD

metawrap annotate_bins -t $cores -o $1 -b $2
