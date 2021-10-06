#!/bin/bash

########## 1 INITIAL PROKARYOTIC BINNING  ###################

# arguments declaration    
assembly=$1                     # assembly fasta file
forward_library=$2              # forward library path
reverse_library=$3              # reverse library path
output_folder=$4                # output folder to be created inside master output folder
num_cores=$5                    # number of threads
memory=$6

echo 'START mudoger module 4-2'

echo 'assembly',assembly
echo 'library 1:', forward_library
echo 'library 2:', reverse_library
echo 'output folder:',output_folder
echo 'num cores:',num_cores

big_bins="$output_folder"/eukrep/filtered_bin_size
annotation="$output_folder"/annotation

echo big_bins, $big_bins

mkdir -p "$annotation"

for bin in "$big_bins"/*; 
do echo $bin; 
mkdir -p "$annotation"/"$(echo "$bin" | rev | cut -f1 -d'/' | rev  | sed "s/.fa//g")";
cp $bin "$annotation"/"$(echo "$bin" | rev | cut -f1 -d'/' | rev  | sed "s/.fa//g")";
cd "$annotation"/"$(echo "$bin" | rev | cut -f1 -d'/' | rev  | sed "s/.fa//g")";
cd -;
done


echo 'END mudoger module 4-2'
