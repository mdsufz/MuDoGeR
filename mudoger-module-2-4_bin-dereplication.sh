#!/bin/bash

# 4 BIN REDUNANCY REMOVAL
cd "$1"

lib=$libname;
bin_count=0; 
mkdir -p unique_bins; 
md5sum refinement-bac/metawrap*bins/*fa  >  unique_bins/md5_sum ;
md5sum refinement-arc/metawrap*bins/*fa  >> unique_bins/md5_sum ;
cat unique_bins/md5_sum | cut -f1 -d' ' | sort | uniq > unique_bins/md5_unique; 

while read l; 
do 
bininit="$(grep "$l" unique_bins/md5_sum  | head -n1 | cut -f3 -d' ')";
binafter=unique_bins/"$lib"-bin."$bin_count".fa; 
cp $bininit $binafter ;
bin_count=$[$bin_count+1];
done < unique_bins/md5_unique ; 
