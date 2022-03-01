#!/bin/bash

# inputs
# 1 is the viral multi fasta file
# 2 is the bins folder
# 3 is the output folder
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database


conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/wish_env


uvigs_file="$1"
bins_folder="$2"
output_folder="$3"

#commands
#1 create output folders
mkdir -p "$output_folder"
mkdir -p "$output_folder"/potential_host_genomes
mkdir -p "$output_folder"/uvigs
mkdir -p "$output_folder"/output_results
mkdir -p "$output_folder"/nullmodels

if [ ! -f "$output_folder"/output_results/prediction.list ];
then
#2 cp and prepare data
yes | cp "$bins_folder"/*fa "$output_folder"/potential_host_genomes
python3 "$(echo $PATH | cut -f1 -d':' | sed "s/conda\/envs\/wish_env\/bin//g" )"/split-all-seq.py "$uvigs_file" "$output_folder"/uvigs/uvig

#3 build viral model
WIsH -c build -g "$output_folder"/potential_host_genomes/ -m "$output_folder"/modelDir

#4 build bacterial null model
WIsH -c predict -g $DATABASES_LOCATION/wish -m "$output_folder"/modelDir -r "$output_folder"/nullmodels -b "$output_folder"/nullmodels

#5 run R script to get nullparameters.tsv
yes | cp "$MUDOGER_DEPENDENCIES_PATH"/computeNullParameters.R  "$output_folder"/nullmodels
cd "$output_folder"/nullmodels
Rscript computeNullParameters.R
cd -

#6 run prediction with null model
WIsH -t 20 -c predict -g "$output_folder"/uvigs/ -m  "$output_folder"/modelDir -r "$output_folder"/output_results/ -b 1 -n "$output_folder"/nullmodels/nullParameters.tsv
else
echo "-> Host prediction is finished"
fi


