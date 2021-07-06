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
forward_library = $3              # forward library path /path/to/libname_1.fastq
reverse_library = $4             # reverse library path. /path/to/libname_2.fastq
output_folder = $5              # master output folder to be created and defined by user
cores=1
