#!/bin/bash


cores="$2"
input_bins_folder="$1"/binning/unique_bins
output_results="$1"/metrics/checkm_qc

# loading conda environment
conda activate metawrap-env

source Mudoger/installation/config_db

extension="fa"
checkm data setRoot "$CHECKM_DB"

# necessary path variable
checkm lineage_wf -t $cores --reduced_tree --tab_table -x $extension -f "$1"/metrics/checkm_qc/outputcheckm.tsv $input_bins_folder $output_results


conda deactivate
