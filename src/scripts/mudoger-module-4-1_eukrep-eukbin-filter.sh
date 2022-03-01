#!/bin/bash

########## 1 INITIAL PROKARYOTIC BINNING  ###################

# fetch eukaryotic contigs/scaffolds
# loading conda environment
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/eukrep_env


# arguments declaration    
assembly=$1                     # assembly fasta file
forward_library=$2              # forward library path
reverse_library=$3              # reverse library path
output_folder=$4                # output folder to be created inside master output folder
num_cores=$5                    # number of threads
memory=$6



EukRep -i "$assembly" --prokarya "$output_folder"/prokaryotic_contigs.fa -o "$output_folder"/eukaryotic_contigs.fa
conda deactivate


# run concoct binnning
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/metawrap_env
metawrap binning -o "$output_folder" -t "$num_cores" -a "$output_folder"/eukaryotic_contigs.fa --concoct "$forward_library" "$reverse_library" # -m "$memory"


# keep bins bigger than 2,5mb
mkdir "$output_folder"/filtered_bin_size
for bin in "$output_folder"/concoct_bins/*fa;
do  size="$(stat $bin | grep Size | cut -f1  | cut -f4 -d' ' )";
if [ $size -gt 2500000 ]; then  cp "$bin" "$output_folder"/filtered_bin_size ;
else : ; fi ;
done





