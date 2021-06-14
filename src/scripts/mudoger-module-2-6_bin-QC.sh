#!/bin/bash


extension="fa"
cores="$1"
input_bins_folder="$2"/binning/unique_bins
output_results="$2"/metrics/checkm_qc

# loading conda environment
conda activate metawrap-env


# necessary path variable
checkm lineage_wf -t $cores  --reduced_tree --tab_table -x $extension -f "$1"/metrics/checkm_qc/outputcheckm.tsv $input_bins_folder $output_results


conda deactivate
