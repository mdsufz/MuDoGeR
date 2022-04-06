#!/bin/bash

#gOTUpick: a tool to separate genomes below species level

#The user should have initial data to use BRAT: output from GTDB-Tk, output from CheckM, output from BBMap, the recovered bins/FASTA in fasta format, and Quality Controlled Pair End reads.

conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database


bbtools_results_path="$1"
checkm_results_path="$2"
gtdbtk_results_path="$3"
mags_results_path="$4"
output_path="$5"
cores="$6"


conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env

mkdir -p $output_path


bash -i gOTUpick.sh --fastANI-thread $cores --bb-input path/to/BBMap-input --checkm-input path/to/CheckM-input --gtdb-input path/to/gtdb-input -m path/to/mags -o path/to/outputdir

