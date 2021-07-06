# help message


forward_library=$1              # forward library path /path/to/libname_1.fastq
reverse_library=$2              # reverse library path. /path/to/libname_2.fastq
output_folder=$3                # master output folder to be created and defined by user
num_cores=$4                    # number of cores to use

lib="$(echo "$forward_library" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | cut -f1 -d'_' )"          # create output master
mkdir -p $lib

master_output_dir="$output_folder"/"$lib"

# 1 QUALITY CONTROL (QC) OF READS
bash -i MuDoGeR/src/scripts/mudoger-module-1-1_QC.sh "$forward_library" "$reverse_library" "$master_output_dir"/qc "$num_cores"

# 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY
bash -i MuDoGeR/src/scripts/mudoger-module-1-2_kmermempred.sh "$master_output_dir"/qc/final_pure_reads_1.fastq "$master_output_dir"/khmer

# 3 ASSEMBLY
bash -i MuDoGeR/src/scripts/mudoger-module-1-3_assembly.sh "$master_output_dir"/qc/final_pure_reads_1.fastq "$master_output_dir"/qc/final_pure_reads_2.fastq  "$master_output_dir" "$num_cores"
