#!/bin/bash

########## RUNNING GENEMARK ###################

# loading conda environment and adjusting the paths
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

#Set genemark scripts folder as variable
genemarker_scripts_folder="$MUDOGER_DEPENDENCIES_ENVS_PATH"/genemarker_env/GENEMARK_MAIN/gmes_linux_64*

# arguments declaration    
euk_folder=$1 #"$libname_folder"/eukaryotes

filtered_euk_bins_folder="$1/filtered_euk_bins"     #folder containing the filtered eukaryotic bins
output_folder="$1/genemarker_annotation"            #output folder


# Run Genemark

#if [ -f  "$output_folder"/eukaryotic_contigs.fa ]; 		# if one of the outputs is already there, do not run
#then echo "-> EukRep already done. Please check here: "$output_folder"/eukaryotic_contigs.fa"
#else

echo -e "\n --->RUNNING GENEMARK"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/genemarker_env
mkdir -p $output_folder

for bin_path in $filtered_euk_bins_folder/*; do

bin=`echo ${bin_path} | rev | cut -f1 -d'/' | rev`

if [ -f  $output_folder/"$bin"_genemark/genemark.gtf ];
then
:
else

mkdir -p $output_folder/"$bin"_genemark
cd $output_folder/"$bin"_genemark
cp $bin_path $output_folder/"$bin"_genemark

perl $genemarker_scripts_folder/gmes_petap.pl  --ES -min_contig 3000 --sequence $output_folder/"$bin"_genemark/"$bin" 
cd -
fi
done
conda deactivate

echo -e "\n --->END GENEMARK"

#End Genemark
