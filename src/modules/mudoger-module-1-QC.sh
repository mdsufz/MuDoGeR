
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
