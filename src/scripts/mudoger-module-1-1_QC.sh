#!/bin/bash

################# STARTING MODULE 1-1  ###############################################
echo '---------- STARTING MODULE 1-1 QC'

# 1 QUALITY CONTROL (QC) OF READS
# 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY
# 3 ASSEMBLY

########## 1 QUALITY CONTROL (QC) OF READS  ###################
## load conda metawrap conda quality control
conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/metawrap_env

# arguments declaration
log="log_qc"                      # definition of path to QC log       
forward_library=$1              # forward library path
reverse_library=$2              # reverse library path
qc_folder=$3                  # output folder to be created inside master output folder
bm_tag="--skip-bmtagger"          # define if human reads should be filtered out or not
num_cores=$4                     # number of threads

# command of quality control
metawrap read_qc "$bm_tag" -1 "$forward_library" -2 "$reverse_library" -t "$num_cores" -o "$qc_folder"  &> "$out_qc"/"$log"

# leave conda environment
conda deactivate
