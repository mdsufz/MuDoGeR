#!/bin/bash

#gOTUpick: a tool to separate genomes below species level

#The user should have initial data to use BRAT: output from GTDB-Tk, output from CheckM, output from BBMap, the recovered bins/FASTA in fasta format, and Quality Controlled Pair End reads.

conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

#Input parameters
WORKDIR="$1"
metadata_table="$2"
cores="$3"

#Define dependent parameters
bbtools_input_path="prokaryotes/metrics/genome_statistics/bbtools.tsv"
checkm_input_path="prokaryotes/metrics/checkm_qc/outputcheckm.tsv"
gtdbtk_input_path="prokaryotes/metrics/GTDBtk_taxonomy/gtdbtk.bac120.summary.tsv"
bins_input_path="prokaryotes/binning/unique_bins/"

#Define output files path
all_bins_path="$WORKDIR/mapping_results/all_bins"
all_metrics_path="$WORKDIR/mapping_results/all_metrics"
gOTUpick_results_path="$WORKDIR/mapping_results/gOTUpick_results"

#gOTU pick started
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env

mkdir -p $all_bins_path
mkdir -p $all_metrics_path
mkdir -p $gOTUpick_results_path

#Concat results metrics for all samples
awk '
    FNR==1 && NR!=1 {next;}{print}
' $WORKDIR/*/$bbtools_input_path >$all_metrics_path/bbtools_all.tsv


awk '
    FNR==1 && NR!=1 {next;}{print}
' $WORKDIR/*/$checkm_input_path >$all_metrics_path/checkm_all.tsv


awk '
    FNR==1 && NR!=1 {next;}{print}
' $WORKDIR/*/$gtdbtk_input_path >$all_metrics_path/gtdbtk_all.tsv


#copy unique bins
yes | cp $WORKDIR/*/$bins_input_path/* $all_bins_path


#Run

bash -i "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env/bin/gOTUpick.sh --fastANI-thread $cores --bb-input $all_metrics_path/bbtools_all.tsv --checkm-input $all_metrics_path/checkm_all.tsv --gtdb-input $all_metrics_path/gtdbtk_all.tsv -m $all_bins_path -o $gOTUpick_results_path

conda deactivate
