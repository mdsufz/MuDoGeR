#!/bin/bash

# loading conda environment
echo '------- START MODULE 2-5 BIN TAXONOMY'
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

# load gtdbtk env
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/gtdbtk_env

GTDBTK_DATA_PATH=$(realpath "$DATABASES_LOCATION"gtdbtk/release*)

export GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH

master_folder=$1
cores=$2

#Run
if [ -f "$master_folder"/metrics/GTDBtk_taxonomy/gtdbtk.bac120.summary.tsv ] || [ -f "$master_folder"/metrics/GTDBtk_taxonomy/gtdbtk.ar53.summary.tsv ];
then echo "";
else
gtdbtk classify_wf --extension fa --cpus "$cores" --genome_dir "$master_folder"/binning/unique_bins --out_dir "$master_folder"/metrics/GTDBtk_taxonomy
fi

#Create merged results file
if [ -f "$master_folder"/metrics/GTDBtk_taxonomy/gtdbtk.bac120.summary.tsv ] || [ -f "$master_folder"/metrics/GTDBtk_taxonomy/gtdbtk.ar53.summary.tsv ]; 
then awk 'FNR==1 && NR!=1 {next;}{print}' "$master_folder"/metrics/GTDBtk_taxonomy/gtdbtk.*summ*.tsv > "$master_folder"/metrics/GTDBtk_taxonomy/gtdbtk_result.tsv;
echo "GTDBtk results generated!"
else
echo "Error: GTDBtk summary files not found"
fi

conda deactivate
