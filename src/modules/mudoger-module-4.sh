#!/bin/bash
# EUKARYOTIC MODULES FILESYSTEM STRUCTURE
# EUKARYA SORTING AND BINNING
#     EUKREP                         submodule 3-1
#     METAWARAP CONCOCT              submodule 3-1
#     SIZE FILTERING                 submodule 3-1
#     GENEMARK (GENE PREDICTION)    submodule 3-2
#     MAKER (GENE ANNOTATION)        submodule 3-3
#     BUSCO (COMPLETENSS COMPUTATION)  submodule 3-4
#     EUKCC (QUALITY ASSESSMENT)
#     
#     METRICS 
#       N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more... submodule 3-5


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


config_path="$(which config.sh)"
source $config_path
# 1 EUKREP CONTIG FILTERING AND  BINNING USING CONCOCT

mkdir -p "$libname_folder"/eukaryotes/eukaryotic_contigs
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-4-1_eukrep-eukbin-filter.sh "$assembly"       \
                                      "$forward_library"                                 \
                                      "$reverse_library"                                 \
                                      "$libname_folder"/eukaryotes			 \
                                      "$cores"                                           \  
				      "$memory"

bash -i MuDoGeR/src/scripts/mudoger-module-4-2_genemark-maker-busco.sh "$assembly"       \
                                      "$forward_library"                                 \
                                      "$reverse_library"                                 \
                                      "$libname_folder"/eukaryotes			 \
                                      "$cores"                                           \  
				      "$memory"



