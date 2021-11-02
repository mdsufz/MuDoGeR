#!/bin/bash


# loading conda environment
conda activate gtdbtk-env

source MuDoGeR/installation/config_db

GTDBTK_DATA_PATH=/mnt/tools/miniconda2/envs/gtdbtk-env/database/release202

master_folder=$1
cores=$2

# necessary path variable
#export GTDBTK_DATA_PATH="$gtdbtk_database"

gtdbtk classify_wf --extension fa --cpus "$cores" --genome_dir "$master_folder"/binning/unique_bins --out_dir "$master_folder"/metrics/GTDBtk_taxonomy

conda deactivate
