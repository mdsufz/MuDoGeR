#!/bin/bash
# PROKARYOTES MODULES FILESYSTEM STRUCTURE
#   BINNING
#     INITIAL-BINNING
#     REFINEMENT-BACTERIA
#     REFINEMENT-ARCHAEA
#     BIN-FINAL-SET (REDUNDANCY REMOVAL)
#   METRICS
#     CHECKM
#     GTDBtk
#     PROKKA
#     STATS (GENOMIC)
#       N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more...



libname_folder=$1              # output folder path 
assembly=$2                    # path to the assembly file
forward_library=$3             # forward library path /path/to/libname_1.fastq
reverse_library=$4             # reverse library path. /path/to/libname_2.fastq
#output_folder=$5              # master output folder to be created and defined by user
cores=$6


#lib="$( echo "$output_folder"/"$(echo "$forward_library" | rev | cut -f2 -d'/' | rev | cut -f1 -d'.' | cut -f1 -d'_' )")"          # create output master



# 1 INITIAL BINNING USING METAWRAP (CONCOCT, METABAT2, MAXBIN2)
mkdir -p "$libname_folder"/prokaryotes/binning
bash -i MuDoGeR/src/scripts/mudoger-module-2-1_initial-binning.sh "$assembly"          \
                                      "$forward_library"   \
                                      "$reverse_library"   \
                                      "$libname_folder"/prokaryotes/binning/initial-binning \
                                      "$cores"


# 2 BIN REFINEMENT USING METAWRAP FOR BACTERIA (50,10)
bash -i MuDoGeR/src/scripts/mudoger-module-2-2_bin-ref-bacteria.sh "$libname_folder"/prokaryotes/binning/refinement-bac                \
                                       "$cores"                                                     \
                                       "$assembly"                                                  \
                                       "$libname_folder"/prokaryotes/binning/initial-binning/concoct_bins  \
                                       "$libname_folder"/prokaryotes/binning/initial-binning/maxbin2_bins  \
                                       "$libname_folder"/prokaryotes/binning/initial-binning/metabat2_bins
              
              

# 3 BIN REFINEMENT USING METAWRAP FOR ARCHAEA (40,30)
bash -i MuDoGeR/src/scripts/mudoger-module-2-3_bin-ref-archea.sh "$libname_folder"/prokaryotes/binning/refinement-arc                \
                                     "$cores"                                                     \
                                     "$assembly"                                                  \
                                     "$libname_folder"/prokaryotes/binning/initial-binning/concoct_bins  \
                                     "$libname_folder"/prokaryotes/binning/initial-binning/maxbin2_bins  \
                                     "$libname_folder"/prokaryotes/binning/initial-binning/metabat2_bins
              
# 4 BIN REDUNANCY REMOVAL
mkdir -p "$libname_folder"/prokaryotes/binning/unique_bins
cd "$libname_folder"/prokaryotes/binning
bash -i MuDoGeR/src/scripts/mudoger-module-2-4_bin-dereplication.sh  "$libname_folder" 
cd -

# 5 GTDBtk taxonomy assignment
mkdir -p "$libname_folder"/prokaryotes/metrics
bash -i MuDoGeR/src/scripts/mudoger-module-2-5_bin-taxonomy.sh "$libname_folder"/prokaryotes "$cores"

# 6 CheckM quality control
bash -i MuDoGeR/src/scripts/mudoger-module-2-6_bin-QC.sh  "$libname_folder"/prokaryotes/ "$cores"

# 7 PROKKA Annotation
bash -i MuDoGeR/src/scripts/mudoger-module-2-7_prokka.sh  "$libname_folder"/prokaryotes/ "$cores"

# 8 Metrics: N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more...
bash -i MuDoGeR/src/scripts/mudoger-module-2-8_genomics-metrics.sh  "$libname_folder"/prokaryotes/ "$cores"

# 9 WRAPUP INTO NICE TABLE

#
