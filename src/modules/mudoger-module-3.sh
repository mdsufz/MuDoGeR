#!/bin/bash
# VIRUSES MODULES FILESYSTEM STRUCTURE
#   VIRAL INVESTIGATION
#     VIRSORTER               submodule 3-1
#     VIRFINDER               submodule 3-1
#     VIBRANT                 submodule 3-1
#     DEREPLICATION           submodule 3-1
#   METRICS
#     VCONTACT2 (TAXONOMY)     submodule 3-2
#     WISH (HOST PREDICTION)   submodule 3-3
#     VCHECK (QUALITY)         submodule 3-4
#     STATS (GENOMIC)
#       N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more... submodule 3-5



#libname=$1
assembly=$1
#forward_library = $3            # forward library path /path/to/libname_1.fastq    # commented because reads are not used
#reverse_library = $4            # reverse library path. /path/to/libname_2.fastq   # and thus not necessary
libname_folder=$2              # master output folder to be created and defined by user
cores=10

# 1 viral investigation
mkdir -p "$libname_folder"/viruses/investigation
bash -i MuDoGeR/src/scripts/mudoger-module-3-1-viral_investigation.sh "$assembly" "$libname_folder"/viruses/investigation $cores

# 2 viral vcontact2 (taxonomy)
bash -i MuDoGeR/src/scripts/mudoger-module-3-2_vcontact2.sh "$libname_folder"/viruses/investigation/dereplication/uvigs_95-70.fna "$libname_folder"/viruses/taxonomy "$cores"

# 3 host prediction
#commented for testing 
#bash -i MuDoGeR/src/scripts/mudoger-module-3-3_host-prediction.sh "$libname_folder"/viruses/investigation/dereplication/uvigs_95-70.fna  "$libname_folder"/prokaryotes/binning/unique_bins "$libname"/viruses/dereplication/host_prediction

# 4 vcheck
#commented for testing
#bash -i MuDoGeR/src/scripts/mudoger-module-3-4_vcheck.sh "$libname_folder"/viruses/dereplication/uvigs/uvigs_95-70.fna "$libname_folder"/viruses/vcheck_quality

# 5 uvigs metrics
# bash -i MuDoGeR/src/scripts/mudoger-module-3-5_uvigs-metrics.sh


