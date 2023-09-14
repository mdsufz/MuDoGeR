#!/bin/bash

conda activate mudoger_env

help_message () {
	echo ""
	echo "Usage: bash -i mudoger-module-1.sh -1 reads_1.fastq -2 reads_2.fastq -o output_dir -t num_cores"
	echo "Options:"
	echo ""
	echo "	-1 STR			forward fastq reads path"
	echo "	-2 STR			reverse fastq reads path"
	echo "	-o STR			output directory path"
	echo "	-m INT			given Memory to the Assembly process in GB (default=10)"
	echo "	-t INT			number of threads/cores (default=1)"
	echo "	--megahit		assemble with megahit (default)"
	echo "	--metaspades		assemble with metaspades instead of megahit"
 	echo "	--m-spec		Flag to check if the user has defined RAM to assemble"
	echo "	--skip-bmtagger		Do not run bmtagger to remove human reads"
	echo "	-h --help		print this message"
	echo "";}

#DEFINE DEFAULT PARAMETERS

# The option memory was disabled as we are using the memory predicted by script 1.2. We need to develop another way if it does not work with the amount of memory predicted. 

output_folder=$(pwd) 		#output path for the downloaded sequences
num_cores=1 			#number of threads that are going to be used
megahit=""			#assemble with megahit (default)"
metaspades=""			#assemble with metaspades instead of megahit"
bm_tag=""			#Do not run bmtagger


# loop through input params
while true; do
	case "$1" in
		-1) forward_library=$2; shift 2;;
		-2) reverse_library=$2; shift 2;;
		-o) output_folder=$2; shift 2;;
		-t) num_cores=$2; shift 2;;
		-m) memory=$2; shift 2;;
		--metaspades) metaspades="--metaspades"; shift 1;;
  		--m-spec) memory_specified=$2; shift 2;;
		--skip-bmtagger) bm_tag="--skip-bmtagger"; shift 1;;
		-h | --help) help_message; exit 1; shift 1;;
		--) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done


master_output_dir="$output_folder"

mkdir -p "$master_output_dir"

config_path="$(which config.sh)"
source $config_path


# 1 QUALITY CONTROL (QC) OF READS
if [ -f  "$master_output_dir"/qc/final_pure_reads_1.fastq ]; 		# if one of the outputs is already there, do not run
then echo "-> Reads were already qc'ed. Please check here: "$master_output_dir"/qc"
else echo "-> Running QC"
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-1-1_QC.sh "$forward_library" "$reverse_library" "$master_output_dir"/qc "$num_cores" "$bm_tag"
fi

# 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY
if [ "$memory_specified" = "true" ]; then
    echo "User has specified memory as $memory GB. Skipping memory prediction."
else
    echo "User has not specified memory."
    if [ -f "$master_output_dir/khmer/final_prediction.tsv" ]; then
        echo "-> Memory prediction is done. Please check here: $master_output_dir/khmer"
    else
        echo "-> Running Khmer pred"
        bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-1-2_kmermempred.sh "$master_output_dir"/qc/final_pure_reads_1.fastq "$master_output_dir"/khmer
    fi
fi

# 3 ASSEMBLY
if [ -f "$master_output_dir/assembly/final_assembly.fasta" ]; then
    echo "-> Assembly is done. Please check here: $master_output_dir/assembly"
else
    echo "-> Running assembly"
    
    if [ "$memory_specified" = "true" ]; then
        echo "User has specified memory as $memory GB. Using user-specified memory."
    else
        if [ -f "$master_output_dir/khmer/final_prediction.tsv" ]; then
            echo "Reading memory requirements from $master_output_dir/khmer/final_prediction.tsv"
            mem_mb="$(tail -n1 "$master_output_dir/khmer/final_prediction.tsv" | cut -f2)"
            mem_gb="$((mem_mb / 1000))"
            memory=$mem_gb
            echo "Using predicted memory for assembly: $memory GB"
        else
            echo "Memory not specified and prediction file not found."
            echo "Please either specify memory with -m or make sure the previous step was done successfully."
            echo "Using default memory of 100 GB."
            memory=100
        fi
    fi

bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-1-3_assembly.sh "$master_output_dir"/qc/final_pure_reads_1.fastq "$master_output_dir"/qc/final_pure_reads_2.fastq  "$master_output_dir" "$num_cores" "$metaspades" "$memory"
fi

echo "-> module 1 (preprocess and assembly) is finished"

conda deactivate
