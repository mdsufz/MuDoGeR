# 1 QUALITY CONTROL (QC) OF READS
# 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY
# 3 ASSEMBLY


################# STARTING MODULE 1  ###############################################
output_master="output_folder"
mkdir -p "$output_master"       # create master output folder

########## 1 QUALITY CONTROL (QC) OF READS  ###################
## load conda metawrap conda quality control
conda activate metawrap-env

# arguments declaration
log="log_qc"                      # definition of path to QC log       
forward_library = $1              # forward library path
reverse_library = $2              # reverse library path
output_folder = $3                # output folder to be created inside master output folder
bm_tag="--skip-bmtagger"          # define if human reads should be filtered out or not
num_cores = 1                     # number of threads

# command of quality control
metawrap read_qc "$bm_tag" -1 "$forward_library" -2 "$reverse_library" -t "$num_cores" -o 
"$output_master"/"$output_folder"  &> "$log"

# leave conda environment
conda deactivate


################### 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY ###################
## load conda metawrap conda quality control
conda activate conda-env  # conda to be created
# https://khmer.readthedocs.io/en/v2.1.1/user/install.html

# arguments declaration
lib="$output_master"/"$output_folder"/final_pure_reads_1.fastq # fastq file to be investigated after qc
kmer_mem_pred = "$output_master"/kmer_mem_predict              # folder where the output will be dumped
mkdir -p "$kmer_mem_pred"

kmers_script="/path/to/unique-kmers.py"    # to be defined where

# run khmer for sizes 33 and 55
python3 "$kmers_script" -k 33 -R "$kmer_mem_pred"/kmer-33 "$lib"
python3 "$kmers_script" -k 55 -R "$kmer_mem_pred"/kmer-55 "$lib"

# extract results from khmer and dump into output file
echo -e "$(echo $1 | rev | cut -f1 -d'/' | rev )\t\c" > "$kmer_mem_pred"/input_for_predictR.tsv
echo -e "$(cat "$2"/kmer-33 | head -n1 | cut -f1 -d' ')\t\c" >> "$kmer_mem_pred"/input_for_predictR.tsv
echo -e "$(cat "$2"/kmer-55 | head -n1 | cut -f1 -d' ')" >> "$kmer_mem_pred"/input_for_predictR.tsv

# copy necessary files for memory prediction. those scripts in-house and should be downloadable via git
cp path/to/function_predict_memory.R  "$kmer_mem_pred"
cp path/to/models.RData  "$kmer_mem_pred"
cp path/to/predict.R  "$kmer_mem_pred"

# run memory prediction
cd  "$kmer_mem_pred"                                                # move to the kmer output folder
Rscript predict.R input_for_predictR.tsv                            # run memory prediction
cat metaspades_prediction.tsv | cut -f1,10 > final_prediction.tsv   # parse to final file
cd -                                                                # move back to previous folder

# erase auxiliary files
rm "$kmer_mem_pred"/models.RData
rm "$kmer_mem_pred"/function_predict_memory.R
rm "$kmer_mem_pred"/predict.R

# finish
echo “end predict_mem_input”
date
# NOTE: the outputted memory prediction is in mega bytes. 
# if the assembly resources request is done using gigabytes
# conversion of megabytes to gigabytes
mem_mb="$(tail -n1 "$kmer_mem_pred"/final_prediction.tsv    | cut -f2 | cut -f2 -d' ')"
mem_gb="$(echo $((mem_mb / 1000)))"


################### 3 ASSEMBLY  #########################################################
## load conda metawrap conda quality control
conda activate metawrap-env

assembly_folder="$output_master"/ASSEMBLY
mkdir -p "$assembly_folder"      # create assembly output folder

# arguments declaration
log="log_assembly"                # definition of path to assembly log       
forward_library = $1              # forward library path
reverse_library = $2              # reverse library path
output_folder = $3                # output folder to be created inside master output folder
num_cores = 1                     # number of threads
assembler="--metaspades"          # --metaspades or --megahit
# assembly command
metawrap assembly -1 $forward_library -2 $reverse_library -m $mem_gb -t "$num_cores" "$assembler" -o "$assembly_folder"/"$output_folder"

if [ -f "$assembly_folder"/"$output_folder"/final_assembly.fa ]; 
then echo  "assembly was succesful" ; 
else echo "assembly didnt work. please check your resources"; 
fi


