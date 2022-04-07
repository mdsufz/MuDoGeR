#!/bin/bash

########## RUNNING EUKCC ###################

# loading conda environment and adjusting the paths
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

#Set eukcc database location
eukcc_db_path="$DATABASES_LOCATION/eukccdb"


# arguments declaration    
euk_folder=$1 #"$libname_folder"/eukaryotes

filtered_euk_bins_folder="$1/filtered_euk_bins"     #folder containing the filtered eukaryotic bins
output_folder="$1/eukcc_quality"            #output folder

num_core=$2


# Run EUKCC

echo -e "\n --->RUNNING EUKCC"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/eukcc_env
mkdir -p $output_folder

for bin_path in $filtered_euk_bins_folder/*; do

bin=`echo ${bin_path} | rev | cut -f1 -d'/' | rev`

if [ -f  $output_folder/"$bin"_eukcc/eukcc.csv ];
then
:
else

mkdir -p $output_folder/"$bin"_eukcc


eukcc single --db $eukcc_db_path --threads $num_core --out $output_folder/"$bin"_eukcc $bin_path

fi
done
conda deactivate

echo -e "\n --->END EUKCC QUALITY CALCULATION"

#End EUKCC
