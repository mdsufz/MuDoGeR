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
prefilter="--no-prefilter"

#Define dependent parameters
bbtools_input_path="prokaryotes/metrics/genome_statistics/bbtools.tsv"
checkm_input_path="prokaryotes/metrics/checkm_qc/outputcheckm.tsv"
gtdbtk_input_path="prokaryotes/metrics/GTDBtk_taxonomy/gtdbtk_result.tsv"
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

awk -F$'\t' '{print $NF}' $all_metrics_path/bbtools_all.tsv | rev | cut -d '/' -f1 | sed 's/af.//1' | rev > $all_metrics_path/aux

awk -F"\t" '{OFS=FS}{ $20="" ; print }' $all_metrics_path/bbtools_all.tsv > $all_metrics_path/bbtools_all_tmp.tsv

paste --delimiters='\t' $all_metrics_path/bbtools_all_tmp.tsv $all_metrics_path/aux > $all_metrics_path/bbtools_all_final.tsv

mv -f $all_metrics_path/bbtools_all_final.tsv $all_metrics_path/bbtools_all.tsv

rm -f $all_metrics_path/bbtools_all_tmp.tsv
rm -f $all_metrics_path/aux
rm -f $all_metrics_path/bbtools_all_final.tsv

awk '
    FNR==1 && NR!=1 {next;}{print}
' $WORKDIR/*/$checkm_input_path >$all_metrics_path/checkm_all.tsv


awk '
    FNR==1 && NR!=1 {next;}{print}
' $WORKDIR/*/$gtdbtk_input_path >$all_metrics_path/gtdbtk_all.tsv


#copy unique bins
yes | cp $WORKDIR/*/$bins_input_path/* $all_bins_path


#Run

bash -i "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env/bin/gOTUpick.sh --fastANI-thread $cores $prefilter --bb-input $all_metrics_path/bbtools_all.tsv --checkm-input $all_metrics_path/checkm_all.tsv --gtdb-input $all_metrics_path/gtdbtk_all.tsv -m $all_bins_path -o $gOTUpick_results_path --a2 95

#Create auxiliary results files

if [ "$(ls $gOTUpick_results_path/tax_groups/ | wc -l)" -eq 1 ]; then

conda deactivate

else
awk 'NR==1 {printf("%s\t%s\n", $0, "representative_bin")}  NR>1 {printf("%s\t%s\n", $0, "*") }' $gOTUpick_results_path/final_output/bestbins.txt > $gOTUpick_results_path/final_output/repbin_aux

awk '{ print $2, "\t" ,$4}' $gOTUpick_results_path/results/final_groups.tsv | tail -n +2 > $gOTUpick_results_path/final_output/aux_final_groups

cat $gOTUpick_results_path/final_output/repbin_aux $gOTUpick_results_path/final_output/aux_final_groups > $gOTUpick_results_path/final_output/concat_file.tsv

awk 'BEGIN{OFS=","} {$1=$1; print}' $gOTUpick_results_path/final_output/concat_file.tsv > $gOTUpick_results_path/final_output/concat_file.csv

awk -F "\"*,\"*" '!seen[$1,$2]++' $gOTUpick_results_path/final_output/concat_file.csv | sort -k2 -t, > $gOTUpick_results_path/final_output/final_groups_output.csv

rm -f $gOTUpick_results_path/final_output/repbin_aux
rm -f $gOTUpick_results_path/final_output/aux_final_groups
rm -f $gOTUpick_results_path/final_output/concat_file.tsv
rm -f $gOTUpick_results_path/final_output/concat_file.csv
rm -fr $all_metrics_path

conda deactivate

fi





