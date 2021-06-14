#!/bin/bash


# loading conda environment
conda activate gtdbtk-env


# necessary path variable
source MuDoGeR/installation/config

gtdbtk  classify_wf --extension  fa  --cpus ${NSLOTS:-1} --genome_dir "$1"/binning/unique_bins --out_dir "$1"/GTDBtk_taxonomy


conda deactivate
