#!/bin/bash

########## RUNNING BUSCO ###################

# loading conda environment and adjusting the paths
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

#Set busco database location
busco_db_path="$DATABASES_LOCATION/buscodbs"


# arguments declaration    
euk_folder=$1 #"$libname_folder"/eukaryotes

maker2_genes_folder="$1/maker2_gene_annotation"     #folder containing the euk_bin with any annotated maker2 genes
output_folder="$1/euk_completeness"                 #output folder

num_core=$2


# Run BUSCO

echo -e "\n --->RUNNING BUSCO"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/busco_env
mkdir -p $output_folder
cd $output_folder
for bin_path in $maker2_genes_folder/*; do

bin=`echo ${bin_path} | rev | cut -f1 -d'/' | cut -d '.' -f2- | rev`

if [ -f  $output_folder/"$bin"_busco/run_lineage_name/full_table.tsv ];
then
:
else

#mkdir -p $output_folder/"$bin"_busco/
#cd $output_folder/"$bin"_busco/
sequence_path=$bin_path/"$bin".maker.output/OUTPUT.all.maker.genemark.proteins.fasta

busco -i $sequence_path -l "$busco_db_path"/lineages/eukaryota_odb10 -m protein -o "$bin"_busco -c $num_core

fi
done
cd -
conda deactivate

echo -e "\n --->END BUSCO CALCULATION"

#End BUSCO
