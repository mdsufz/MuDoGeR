#!/bin/bash

# this bash script calculates the reduced BRAT

# loading conda environment
echo '------- START MODULE 5-2 REDUCED BRAT'
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

#Input data
project_folder=$1
meta_table=$2
cores=$3

#Process inputs
output_folder=$project_folder/mapping_results/reduced_BRAT

# load brat env
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/brat_env

mkdir -p $output_folder

#Calculate genome size
mkdir -p $project_folder/mapping_results/genome_size

for bin_path in $project_folder/mapping_results/all_bins/*.fa;
  do  bin="$(echo "$bin_path" | rev | cut -f1 -d'/' | rev )";
  echo -e "\nCalculating genome size from $bin"
  cat "$bin_path" | grep -v ">" > $project_folder/mapping_results/genome_size/aux;
  size="$(wc -c $project_folder/mapping_results/genome_size/aux | cut -f1 -d' ' )";
  echo ${size} ${bin/.fa/} > "$project_folder/mapping_results/genome_size/${bin/.fa/}.genome-size"
  rm -f $project_folder/mapping_results/genome_size/aux;
  
done 

#Indexing bins
cd $project_folder/mapping_results/all_bins/

for bin in *.fa ; 
  do echo -e "\nIndexing $bin" 
  bowtie2-build "$bin" "${bin/.fa/}";
done
cd -




#deactivate env
conda deactivate
