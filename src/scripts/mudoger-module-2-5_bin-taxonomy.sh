#!/bin/bash


# loading conda environment
conda activate gtdbtk-env

source MuDoGeR/installation/config

# necessary path variable
#export GTDBTK_DATA_PATH="$gtdbtk_database"

gtdbtk  classify_wf --extension  fa  --cpus "$2" -genome_dir "$1"/binning/unique_bins --out_dir "$1"/metrics/GTDBtk_taxonomy

conda deactivate
