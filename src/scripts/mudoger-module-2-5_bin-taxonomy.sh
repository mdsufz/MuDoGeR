#!/bin/bash

# loading conda environment
echo '------- START MODULE 2-5 BIN TAXONOMY'
## load conda metawrap conda quality control
conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/gtdbtk-env

source MuDoGeR/installation/config_db

#Run only once during database installation configuration
echo ${CHECKM_DB} | checkm data setRoot ${CHECKM_DB}


GTDBTK_DATA_PATH=/mnt/tools/miniconda2/envs/gtdbtk-env/database/release202

master_folder=$1
cores=$2

# necessary path variable
#export GTDBTK_DATA_PATH="$gtdbtk_database"

gtdbtk classify_wf --extension fa --cpus "$cores" --genome_dir "$master_folder"/binning/unique_bins --out_dir "$master_folder"/metrics/GTDBtk_taxonomy

conda deactivate
