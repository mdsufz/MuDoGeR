#!/bin/bash

# loading conda environment
echo '------- START MODULE 2-6 BIN Quality - CheckM'
# loading conda environment
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database


conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/metawrap_env

#Arguments declaration
output_folder=$1/metrics/checkm_qc/                # output folder to be created inside master output folder
cores=$2
input_bins_folder=$1/binning/unique_bins
extension="fa"

#Create output_folder
mkdir -p $output_folder

#Set CheckM DB
CHECKM_DB="$DATABASES_LOCATION"/checkm

echo ${CHECKM_DB} | checkm data setRoot ${CHECKM_DB}

# necessary path variable
checkm lineage_wf -t $cores --reduced_tree --tab_table -x $extension -f "$output_folder"/outputcheckm.tsv $input_bins_folder $output_folder


conda deactivate
