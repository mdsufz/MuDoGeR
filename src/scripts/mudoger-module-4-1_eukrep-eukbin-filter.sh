
#!/bin/bash

########## 1 INITIAL PROKARYOTIC BINNING  ###################
## load conda metawrap conda quality control
conda activate metawrap-env

# arguments declaration    
assembly=$1                     # assembly fasta file
forward_library=$2              # forward library path
reverse_library=$3              # reverse library path
output_folder=$4                # output folder to be created inside master output folder
num_cores=$5                    # number of threads
memory=$6


# fetch eukaryotic contigs/scaffolds
conda activate eukrep-env
EukRep -i "$assembly" --prokarya "$output_folder"/prokaryotic_contigs -o "$output_folder"/eukaryotic_contigs
conda deactivate


# run concoct binnning
conda activate metawrap-env
metawrap binning -o "$output_folder" -t "$num_cores" -a "$output_folder"/eukaryotic_contigs --concoct "$forward_library" "$reverse_library" # -m "$memory" 


# keep bins bigger than 2,5mb
mkdir filtered_bin_size 
for bin in concoct_binning/concoct_bins/*fa; 
do  size="$(stat $bin | grep Size | cut -f1  | cut -f4 -d' ' )"; 
if [ $size -gt 2500000 ]; then  cp "$bin" filtered_bin_size ; 
else : ; fi ; 
done



