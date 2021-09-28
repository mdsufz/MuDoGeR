#!/bin/bash

# inputs
# 1 is the viral multi fasta file
# 2 is the bins folder
# 3 is the output folder

uvigs_file="$1"
bins_folder="$2"
output_folder="$3"

#commands
#1 create output folders
mkdir -p "$output_folder"   
mkdir -p "$output_folder"/potential_host_genomes
mkdir -p "$output_folder"/viral_particles
mkdir -p "$output_folder"/output_results
mkdir -p "$output_folder"/nullmodels

#2 cp and prepare data
cp "$bins_folder"/*fa "$output_folder"/potential_host_genomes
python MuDoGeR/tools/split-all-seq.py "$uvigs_file" "$output_folder"/viral_particles/viral-particle

#3 build viral model
/data/msb/tools/wish/WIsH/WIsH -c build -g "$output_folder"/potential_host_genomes/ -m "$output_folder"/modelDir

#4 build bacterial null model
/data/msb/tools/wish/WIsH/WIsH -c predict -g /path/to/splittedGenomes_20201218 -m "$output_folder"/modelDir -r "$output_folder"/nullmodels -b     

#5 run R script to get nullparameters.tsv
cp /gpfs1/data/msb/tools/wish/WIsH/computeNullParameters.R  "$output_folder"/nullmodels
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.0
cd "$output_folder"/nullmodels
Rscript computeNullParameters.R 

#6 run prediction with null model
/data/msb/tools/wish/WIsH/WIsH -t 20 -c predict -g "$output_folder"/viral_particles/ -m  "$output_folder"/modelDir -r "$output_folder"/output_results/ -b -n "$output_folder"/nullmodels/nullParameters.tsv


