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



echo 'END mudoger module 4-2'
