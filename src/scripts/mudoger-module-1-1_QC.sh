
# 1 QUALITY CONTROL (QC) OF READS
# 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY
# 3 ASSEMBLY


################# STARTING MODULE 1  ###############################################


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


out_qc="$( echo "$output_folder"/"$(echo "$forward_library" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | cut -f1 -d'_' )")/qc"          # create output master

#mkdir -p $out_qc

# command of quality control
metawrap read_qc "$bm_tag" -1 "$forward_library" -2 "$reverse_library" -t "$num_cores" -o "$out_qc"  &> "$log"

# leave conda environment
conda deactivate
