#!/bin/bash


# download database
# wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz
# conda or virtualenv activation

source /gpfs1/data/msb/tools/GTDB/gtdbtk-v1.3.0/bin/activate

# necessary path variable
export GTDBTK_DATA_PATH=/path/to/release95


gtdbtk  classify_wf --extension  fa  --cpus ${NSLOTS:-1} --genome_dir "$1"/binning/unique_bins --out_dir "$1"/taxonomy


