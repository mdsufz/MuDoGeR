#!/bin/bash

conda activate mudoger_env
help_message () {
	echo ""
	echo "Usage: bash -i MuDoGeR/src/modules/mudoger-module-1.sh -1 reads_1.fastq -2 reads_2.fastq -o output_dir"
	echo "Options:"
	echo ""
	echo "	-1 STR			forward fastq reads path"
	echo "	-2 STR			reverse fastq reads path"
	echo "	-o STR			output directory path"
#	echo "	-m INT			given Memory to the Assembly process in GB (default=10)"
	echo "	-t INT			number of threads/cores (default=1)"
	echo "	--megahit		assemble with megahit (default)"
	echo "	--metaspades		assemble with metaspades instead of megahit"
	echo "	-h --help		print this message"
	echo "";}

#DEFINE DEFAULT PARAMETERS

# the option memory was desabled as we are using the memory predicted by script 1.2. We need to develop another way if it does not work with the amount of memory predicted. 

output_folder=$(pwd) 		#output path for the downloaded sequences
#memory=10			#given Memory to the Assembly process in GB
num_cores=1 			#number of threads that is going to be used
megahit=""			#assemble with megahit (default)"
metaspades="--metaspades"			#assemble with metaspades instead of megahit"


# loop through input params
while true; do
	case "$1" in
		-1) forward_library=$2; shift 2;;
		-2) reverse_library=$2; shift 2;;
		-o) output_folder=$2; shift 2;;
		-t) num_cores=$2; shift 2;;
		#-m) memory=$2; shift 2;;
		--metaspades) metaspades="--metaspades"; shift 1;;
		-h | --help) help_message; exit 1; shift 1;;
		--) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done


#lib="$(echo "$forward_library" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | cut -f1 -d'_' )"          # create output master

master_output_dir="$output_folder"

mkdir -p "$master_output_dir"

config_path="$(which config.sh)"
source $config_path

########## CHANGE the path to the scripts below to a variable Ex: "$mudoger"/src...

# 1 QUALITY CONTROL (QC) OF READS
if [ -f  "$master_output_dir"/qc/final_pure_reads_1.fastq ]; 		# if one of the outputs is already there, do not run
then echo "-> Reads were already qc'ed. Please check here: "$master_output_dir"/qc"
else echo "-> Running QC"
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-1-1_QC.sh "$forward_library" "$reverse_library" "$master_output_dir"/qc "$num_cores"
fi

# 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY
if [ -f  "$master_output_dir"/khmer/metaspades_prediction.tsv ]; 
then echo "-> Memory prediction is done. Please check here: "$master_output_dir"/khmer"
else echo "-> Running khmer pred"
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-1-2_kmermempred.sh "$master_output_dir"/qc/final_pure_reads_1.fastq "$master_output_dir"/khmer
fi



# 3 ASSEMBLY
if [ -f  "$master_output_dir"/assembly/final_assembly.fasta ]; 
then echo "-> Assembly is reading. Please check here: "$master_output_dir"/assembly"
else echo "-> Running assembly"
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-1-3_assembly.sh "$master_output_dir"/qc/final_pure_reads_1.fastq "$master_output_dir"/qc/final_pure_reads_2.fastq  "$master_output_dir" "$num_cores" "$metaspades" #"$memory"
fi

echo "-> module 1 (preprocess and assembly) is finished"

conda deactivate