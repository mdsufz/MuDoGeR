#!/bin/bash

########## RUNNING MAKER2 ###################

# loading conda environment and adjusting the paths
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

#Set genemark scripts folder as variable
maker2_environment="$MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env"

# arguments declaration    
euk_folder=$1 #"$libname_folder"/eukaryotes

genemark_results_folder="$1/genemarker_annotation"     #folder containing the filtered eukaryotic bins
output_folder="$1/maker2_gene_annotation"            #output folder

cores=$2

# Run MAKER2

echo -e "\n --->RUNNING MAKER"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/maker2_env

export ZOE="$maker2_environment/MAKER_MAIN/maker/exe/snap/Zoe"
export AUGUSTUS_CONFIG_PATH="$maker2_environment/config"
export PATH=$maker2_environment/MAKER_MAIN/maker/bin:$PATH

mkdir -p $output_folder

for genemark_bin_result_path in $genemark_results_folder/*; do

bin=`echo "$(basename $genemark_bin_result_path)" | cut -f1 -d '_'`
genemark_gmhm_file=$genemark_bin_result_path/gmhmm.mod
bin_file_path=$genemark_bin_result_path/$bin
#if [ -f  $output_folder/"$bin"_genemark/genemark.gtf ];
#then
#:
#else

mkdir -p $output_folder/"$bin"_maker2
cd $output_folder/"$bin"_maker2
cp $bin_file_path $output_folder/"$bin"_maker2

maker -CTL
sed -i 's|genome=|genome='$bin'|' $output_folder/"$bin"_maker2/maker_opts.ctl #do it only for first occurance
sed -i 's|gmhmm=|gmhmm='$genemark_gmhm_file'|g' $output_folder/"$bin"_maker2/maker_opts.ctl
sed -i 's|model_org=all|model_org=|g' $output_folder/"$bin"_maker2/maker_opts.ctl

maker -g $bin -c $cores

cd -

cd $output_folder/"$bin"_maker2/*output

fasta_merge -d *_master_datastore_index.log -o OUTPUT

cd -
#fi
done
conda deactivate

echo -e "\n --->END MAKER"

#End MAKER
