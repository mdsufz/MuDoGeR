#!/bin/bash

################### 3 ASSEMBLY  #########################################################
## load conda metawrap conda quality control
conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/metawrap_env


# arguments declaration
log="log_assembly"             # definition of path to assembly log       
forward_library=$1             # forward library path
reverse_library=$2             # reverse library path
master_folder=$3               # output folder to be created inside master output folder
num_cores=$4                   # number of threads
assembler=$5                   # --metaspades or --megahit
mem_gb=$6                      # add the amount of memory in case you do not want to use the amount predicted. Give amount in GB. Ex: "10" for 10GB.

assembly_folder="$master_folder"/assembly

mkdir -p "$assembly_folder"      # create assembly output folder


#Define memory requirements
kmer_file="$master_folder"/khmer/final_prediction.tsv

if [ -z "$6" ];
then mem_mb="$(grep "$(echo "$forward_library" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | cut -f1 -d'_')" $kmer_file | cut -f2)"; mem_gb="$(echo $((mem_mb / 1000)))";
else mem_gb="$6"
fi

if [ -z "$mem_gb" ]; 
then mem_gb=100;
fi


# assembly command
metawrap assembly -1 $forward_library -2 $reverse_library -m $mem_gb -t "$num_cores" -o "$assembly_folder" "$assembler"

if [ -f "$assembly_folder"/final_assembly.fasta ]; 
then echo "--> Assembly was successful"; 
else echo "--> Assembly didn't work. Please check your resources"; 
fi
