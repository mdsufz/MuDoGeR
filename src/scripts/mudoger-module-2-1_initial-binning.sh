#!/bin/bash

##########  INITIAL PROKARYOTIC BINNING  ###################
echo '------- START MODULE 2-1 INITIAL BINNING'
## load conda metawrap conda quality control
conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/metawrap_env

# arguments declaration    
assembly=$1                     # assembly fasta file
forward_library=$2              # forward library path
reverse_library=$3              # reverse library path
output_folder=$4                # output folder to be created inside master output folder
num_cores=$5                    # number of threads
memory=$6

#mkdir $prok/bact_ref
#outb="$prok/bact_ref"
#mkdir $prok/arch_ref
#outa="$prok/arch_ref"

#con="$prok/concoct_bins"
#met="$prok/metabat2_bins"
#max="$prok/maxbin2_bins"

metawrap binning -o "$output_folder" -t "$num_cores" -a "$assembly" --metabat2 --maxbin2 --concoct "$forward_library" "$reverse_library" # -m "$memory" 
