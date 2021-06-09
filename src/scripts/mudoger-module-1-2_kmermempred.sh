################### 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY ###################
## load conda metawrap conda quality control
conda activate khmer-env  # conda to be created
# https://khmer.readthedocs.io/en/v2.1.1/user/install.html

out_kmer="$( echo "$output_folder"/"$(echo "$forward_library" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | cut -f1 -d'_' )")/kmer"          # create output master



# arguments declaration
lib="$output_master"/"$output_folder"/final_pure_reads_1.fastq # fastq file to be investigated after qc
kmer_mem_pred = "$out_kmer"              # folder where the output will be dumped
mkdir -p "$kmer_mem_pred"

# run khmer for sizes 33 and 55
python3 unique-kmers.py -k 33 -R "$kmer_mem_pred"/kmer-33 "$lib"
python3 unique-kmers.py -k 55 -R "$kmer_mem_pred"/kmer-55 "$lib"

# extract results from khmer and dump into output file
echo -e "$(echo $1 | rev | cut -f1 -d'/' | rev )\t\c" > "$kmer_mem_pred"/input_for_predictR.tsv
echo -e "$(cat "$2"/kmer-33 | head -n1 | cut -f1 -d' ')\t\c" >> "$kmer_mem_pred"/input_for_predictR.tsv
echo -e "$(cat "$2"/kmer-55 | head -n1 | cut -f1 -d' ')" >> "$kmer_mem_pred"/input_for_predictR.tsv

# copy necessary files for memory prediction. those scripts in-house and should be downloadable via git
cp MuDoGeR/tools/mpred_function_predict_memory.R  "$kmer_mem_pred"
cp MuDoGeR/tools/mpred_models.RData  "$kmer_mem_pred"
cp MuDoGeR/tools/mpred_predict.R  "$kmer_mem_pred"

# run memory prediction
cd  "$kmer_mem_pred"                                                # move to the kmer output folder
Rscript predict.R input_for_predictR.tsv                            # run memory prediction
cat metaspades_prediction.tsv | cut -f1,10 > final_prediction.tsv   # parse to final file
cd -                                                                # move back to previous folder

# erase auxiliary files
rm "$kmer_mem_pred"/models.RData
rm "$kmer_mem_pred"/function_predict_memory.R
rm "$kmer_mem_pred"/predict.R

conda deactivate
# finish
echo “end predict_mem_input”
date
# NOTE: the outputted memory prediction is in mega bytes. 
# if the assembly resources request is done using gigabytes
# conversion of megabytes to gigabytes
mem_mb="$(tail -n1 "$kmer_mem_pred"/final_prediction.tsv    | cut -f2 | cut -f2 -d' ')"
mem_gb="$(echo $((mem_mb / 1000)))"
