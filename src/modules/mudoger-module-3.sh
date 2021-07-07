#!/bin/bash
# VIRUS MODULES FILESYSTEM STRUCTURE
# #!/bin/bash
# PROKARYOTES MODULES FILESYSTEM STRUCTURE
#   VIRAL RECOVERY TOOLS
#     VIRSORTER               submodule 3-1
#     VIRFINDER               submodule 3-1
#     VIBRANT                 submodule 3-1
#     DATA CLEAN DEREPLICATION STEP (STAMPEDE CLUSTERGENOMES) submodule 3-1
#   METRICS
#     VCONTACT2 (TAXONOMY)    submodule 3-2
#     VCHECK (QUALITY)        submodule 3-3
#     WISH (HOST PREDICTION)  submodule 3-4
#     STATS (GENOMIC)
#       N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more... submodule 3-5



libname=$1
assembly = $2
forward_library = $3            # forward library path /path/to/libname_1.fastq
reverse_library = $4            # reverse library path. /path/to/libname_2.fastq
output_folder = $5              # master output folder to be created and defined by user
cores=1

# 1 viral investigation
mkdir -p "$libname"/viruses/investigation
bash -i MuDoGeR/src/scripts/mudoger-module-3-1-viral_investigation.sh "$assembly"  "$forward_library"  "$reverse_library"  "$libname"/viruses/investigation $cores

# 2 viral vcontact2 (taxonomy)
bash -i MuDoGeR/src/scripts/mudoger-module-3-2_vcontact2.sh "$libname"/viruses/dereplication/uvigs/uvigs_95-70.fna "$libname"/viruses/taxonomy "$cores"

# 3 host prediction
bash -i MuDoGeR/src/scripts/mudoger-module-3-3_host-prediction.sh "$libname"/viruses/dereplication/uvigs/uvigs_95-70.fna  "$libname"/prokaryotes/binning/unique_bins

# 4 vcheck
bash -i MuDoGeR/src/scripts/mudoger-module-3-4_vcheck.sh 

# 5 uvigs metrics
bash -i MuDoGeR/src/scripts/mudoger-module-3-5_uvigs-metrics.sh


