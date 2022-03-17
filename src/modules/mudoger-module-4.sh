#!/bin/bash
# EUKARYOTIC MODULES FILESYSTEM STRUCTURE
# EUKARYA SORTING AND BINNING
#     EUKREP                          submodule 4-1
#     METAWARAP CONCOCT               submodule 4-1
#     SIZE FILTERING                  submodule 4-1
#     GENEMARK (GENE PREDICTION)      submodule 4-2
#     EUKCC (QUALITY ASSESSMENT)      submodule 4-3
#     MAKER (GENE ANNOTATION)         submodule 4-4
#     BUSCO (COMPLETENSS COMPUTATION) submodule 4-5
#     
#     METRICS 
#       N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more... submodule 4-6


help_message () {
	echo ""
	echo "Usage: bash -i MuDoGeR/src/modules/mudoger-module-4.sh -1 reads_1.fastq -2 reads_2.fastq -o output_dir -a assembly.fasta -t num_cores"
	echo "Options:"
	echo ""
	echo "	-1 STR			forward fastq reads path"
	echo "	-2 STR			reverse fastq reads path"
	echo "	-a STR			assembly path"
	echo "	-o STR			output directory path"
	echo "	-m INT			given Memory to the Assembly process in GB (default=10)"
	echo "	-t INT			number of threads/cores (default=1)"
	echo "	-h --help		print this message"
	echo "";}

# loop through input params
while true; do
	case "$1" in
		-1) forward_library=$2; shift 2;;
		-2) reverse_library=$2; shift 2;;
		-a) assembly=$2; shift 2;;
		-o) libname_folder=$2; shift 2;;
		-t) cores=$2; shift 2;;
		-m) memory=$2; shift 2;;
		-h | --help) help_message; exit 1; shift 1;;
		--) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done

conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path

#     EUKREP                          submodule 4-1
#     METAWARAP CONCOCT               submodule 4-1
#     SIZE FILTERING                  submodule 4-1

mkdir -p "$libname_folder"/eukaryotes/
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-4-1_eukrep-eukbin-filter.sh "$assembly"       \
                                      "$forward_library"                                 \
                                      "$reverse_library"                                 \
                                      "$libname_folder"/eukaryotes			 \
                                      "$cores"                                           \  
				      "$memory"
echo -e "\n EUK BIN CALCULATION DONE"
exit 0 
#     GENEMARK (GENE PREDICTION)      submodule 4-2

bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-4-2_genemark-maker-busco.sh "$assembly"       \
                                      "$forward_library"                                 \
                                      "$reverse_library"                                 \
                                      "$libname_folder"/eukaryotes			 \
                                      "$cores"                                           \  
				      "$memory"
				      
				      
#     EUKCC (QUALITY ASSESSMENT)      submodule 4-3

bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-4-2_genemark-maker-busco.sh "$assembly"       \
                                      "$forward_library"                                 \
                                      "$reverse_library"                                 \
                                      "$libname_folder"/eukaryotes			 \
                                      "$cores"                                           \  
				      "$memory"

#     MAKER (GENE ANNOTATION)         submodule 4-4

bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-4-2_genemark-maker-busco.sh "$assembly"       \
                                      "$forward_library"                                 \
                                      "$reverse_library"                                 \
                                      "$libname_folder"/eukaryotes			 \
                                      "$cores"                                           \  
				      "$memory"

#     BUSCO (COMPLETENSS COMPUTATION) submodule 4-5

bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-4-2_genemark-maker-busco.sh "$assembly"       \
                                      "$forward_library"                                 \
                                      "$reverse_library"                                 \
                                      "$libname_folder"/eukaryotes			 \
                                      "$cores"                                           \  
				      "$memory"




