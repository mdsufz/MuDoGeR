#!/bin/bash

########## 1 INITIAL EUKARYOTIC BINNING  ###################

# fetch eukaryotic contigs/scaffolds
# loading conda environment
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database


# arguments declaration    
assembly=$1                     # assembly fasta file
forward_library=$2              # forward library path
reverse_library=$3              # reverse library path
output_folder=$4                # output folder to be created inside master output folder
num_cores=$5                    # number of threads
memory=$6

# Run Eukrep

if [ -f  "$output_folder"/eukaryotic_contigs.fa ]; 		# if one of the outputs is already there, do not run
then echo "-> EukRep already done. Please check here: "$output_folder"/eukaryotic_contigs.fa"
else

echo -e "\n --->RUNNING EUKREP"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/eukrep_env

EukRep -i "$assembly" --prokarya "$output_folder"/prokaryotic_contigs.fa -o "$output_folder"/eukaryotic_contigs.fa

conda deactivate
echo -e "\n --->END EUKREP"

fi
#End eukrep

# Run concoct binnning

if [ -f  "$output_folder"/.euk_bin_done ]; 		# if one of the outputs is already there, do not run
then echo "-> Eukariotic bins already done. Please check here: "$output_folder"/eukaryotes_bins"
else

echo -e "\n --->RUNNING CONCOCT BINNING"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/metawrap_env
mkdir -p "$output_folder"/eukaryotes_bins
metawrap binning -o "$output_folder"/eukaryotes_bins -t "$num_cores" -a "$output_folder"/eukaryotic_contigs.fa --concoct "$forward_library" "$reverse_library" 
mv "$output_folder"/eukaryotes_bins/concoct_bins/* "$output_folder"/eukaryotes_bins
rm -fr "$output_folder"/eukaryotes_bins/concoct_bins/
rm -fr "$output_folder"/eukaryotes_bins/work_files/
touch "$output_folder"/.euk_bin_done
conda deactivate
echo -e "\n --->END CONCOCT BINNING"

fi

#End CONCOCT binning

#Start Filtering
if [ -f  "$output_folder"/.euk_filter_done ]; 		# if one of the outputs is already there, do not run
then echo "-> EukRep already done. Please check here: "$output_folder"/filtered_euk_bins"
else

echo -e "\n --->RUNNING FILTER"
# keep bins bigger than 1.5mb
mkdir -p "$output_folder"/filtered_euk_bins
for bin in "$output_folder"/eukaryotes_bins/*fa;
do  size="$(stat $bin | grep Size | cut -f1  | cut -f4 -d' ' )";
if [ $size -gt 1500000 ]; then  cp "$bin" "$output_folder"/filtered_euk_bins ;
else : ; fi ;
done
touch "$output_folder"/.euk_filter_done
echo -e "\n --->END FILTER"

fi

#End filtering





