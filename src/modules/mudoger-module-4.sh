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



libname=$1
assembly=$2
forward_library=$3            # forward library path /path/to/libname_1.fastq
reverse_library=$4            # reverse library path. /path/to/libname_2.fastq
output_folder=$5              # master output folder to be created and defined by user
cores=1

# 1 EUKREP CONTIG FILTERING AND  BINNING USING CONCOCT

mkdir -p "$libname_folder"/eukaryotes/eukaryotic_contigs
#bash -i MuDoGeR/src/scripts/mudoger-module-4-1_eukrep-eukbin-filter.sh "$assembly"       \
#                                      "$forward_library"                                 \
#                                      "$reverse_library"                                 \
#                                      "$libname_folder"/eukaryotes			 \
#                                      "$cores"                                           \  
#				      "$memory"

bash -i MuDoGeR/src/scripts/mudoger-module-4-2_genemark-maker-busco.sh "$assembly"       \
                                      "$forward_library"                                 \
                                      "$reverse_library"                                 \
                                      "$libname_folder"/eukaryotes			 \
                                      "$cores"                                           \  
				      "$memory"



