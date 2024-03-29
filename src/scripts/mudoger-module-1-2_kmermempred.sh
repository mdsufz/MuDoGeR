#!/bin/bash

################### 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY ###################
echo '---------- STARTING MODULE 1-2 KMER'
conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/khmer_env


# arguments declaration
fwd_lib="$1"                   # fastq file to be investigated after qc
khmer_folder="$2"              # folder where the output will be dumped
mkdir -p "$khmer_folder"

# run khmer for sizes 33 and 55
unique-kmers.py -k 33 -R "$khmer_folder"/kmer-33 "$fwd_lib"
unique-kmers.py -k 55 -R "$khmer_folder"/kmer-55 "$fwd_lib"

# extract results from khmer and dump into output file
echo -e "$(echo "$fwd_lib" | rev | cut -f1 -d'/' | rev )\t\c" > "$khmer_folder"/input_for_predictR.tsv
echo -e "$(cat "$khmer_folder"/kmer-33 | head -n1 | cut -f1 -d' ')\t\c" >> "$khmer_folder"/input_for_predictR.tsv
echo -e "$(cat "$khmer_folder"/kmer-55 | head -n1 | cut -f1 -d' ')" >> "$khmer_folder"/input_for_predictR.tsv

# copy necessary files for memory prediction. those scripts in-house and should be downloadable via git
cp "$MUDOGER_DEPENDENCIES_PATH"/mpred_function_predict_memory.R  "$khmer_folder"
cp "$MUDOGER_DEPENDENCIES_PATH"/mpred_models.RData  "$khmer_folder"
cp "$MUDOGER_DEPENDENCIES_PATH"/mpred_predict.R  "$khmer_folder"

# run memory prediction
cd  "$khmer_folder"                                                # move to the kmer output folder
Rscript mpred_predict.R input_for_predictR.tsv                      # run memory prediction
cat metaspades_prediction.tsv | cut -f1,10 > final_prediction.tsv   # parse to final file
cd -                                                                # move back to previous folder

# erase auxiliary files
rm -f "$khmer_folder"/mpred_models.RData
rm -f "$khmer_folder"/mpred_function_predict_memory.R
rm -f "$khmer_folder"/mpred_predict.R

conda deactivate
# finish
echo “end predict_mem_input”
date
# NOTE: the outputted memory prediction is in mega bytes. 
# if the assembly resources request is done using gigabytes
# conversion of megabytes to gigabytes
mem_mb="$(tail -n1 "$khmer_folder"/final_prediction.tsv    | cut -f2 )"
mem_gb="$(echo $((mem_mb / 1000)))"
