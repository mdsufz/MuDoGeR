#!/bin/bash


# input 1 is the fastq file to be investivgated
# input 2 is the folder where the output will be dumped

kmers_script="/gpfs1/data/msb/tools/khmer/scripts/unique-kmers.py"

mkdir -p "$2"

# run khmer for sizes 33 and 55
python3 "$kmers_script" -k 33 -R "$2"/kmer-33 "$1"
python3 "$kmers_script" -k 55 -R "$2"/kmer-55 "$1"

# extract results from khmer and dump into output file
echo -e "$(echo $1 | rev | cut -f1 -d'/' | rev )\t\c" > "$2"/input_for_predictR.tsv
echo -e "$(cat "$2"/kmer-33 | head -n1 | cut -f1 -d' ')\t\c" >> "$2"/input_for_predictR.tsv
echo -e "$(cat "$2"/kmer-55 | head -n1 | cut -f1 -d' ')" >> "$2"/input_for_predictR.tsv

# copy necessary files for memory prediction 
cp /gpfs1/data/msb/PEOPLE/bicalho/Mudoger/ADMIN/scripts/manag_scripts/auxiliary/mem_prediction_for_assembly/function_predict_memory.R /gpfs1/data/msb/PEOPLE/bicalho/Mudoger/ADMIN/scripts/manag_scripts/auxiliary/mem_prediction_for_assembly/models.RData /gpfs1/data/msb/PEOPLE/bicalho/Mudoger/ADMIN/scripts/manag_scripts/auxiliary/mem_prediction_for_assembly/predict.R "$2"

# run memory prediction
cd "$2"
Rscript predict.R input_for_predictR.tsv
cat metaspades_prediction.tsv | cut -f1,10 > final_prediction.tsv
cd -

# erase auxiliary files
rm "$2"/models.RData
rm "$2"/function_predict_memory.R
rm "$2"/predict.R


