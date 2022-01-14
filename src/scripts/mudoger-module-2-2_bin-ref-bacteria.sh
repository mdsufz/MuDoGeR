#!/bin/bash
########## 2 PROKARYOTIC REFINEMENT FOR BACTERIA (50,10  ###################
## load conda metawrap conda quality control
echo '------- START MODULE 2-2 BIN REFINEMENT FOR BACTERIA'
## load conda metawrap conda quality control
conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/metawrap_env



# arguments declaration
#log="log_qc"                      # definition of path to QC log       

output_folder=$1                # output folder to be created inside master output folder
cores=$2
assembly=$3
con=$4
max=$5
met=$6
memory=$7

source MuDoGeR/installation/config_db

#Run only once during database installation configuration
#checkm data setRoot "$CHECKM_DB"

metawrap bin_refinement -o "$output_folder" -t $cores -A "$con" -B "$met" -C "$max" -c 50 -x 10 -m "$memory"

conda deactivate

