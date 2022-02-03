#!/bin/bash


# loading conda environment
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/prokka_env



cores="$2"
input_bins_folder="$1"/binning/unique_bins
output_results="$1"/metrics/prokka


mkdir -p "$output_results"



# for each bin, run prokka and dump results inside output folder
for bin in "$input_bins_folder"/*fa ; do prokka "$bin" --cpus "$cores" --outdir "$output_results"/"$(echo "$bin" | rev | cut -f1 -d'/' | rev | sed "s/.fa//g")"; done

conda deactivate
