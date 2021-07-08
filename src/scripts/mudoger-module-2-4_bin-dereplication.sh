#!/bin/bash

# 4 BIN REDUNANCY REMOVAL

lib=$1;
bin_count=0; 

md5sum "$1"/prokaryotes/binning/refinement-bac/metawrap*bins/*fa  >  "$1"/prokaryotes/binning/unique_bins/md5_sum ;
md5sum "$1"/prokaryotes/binning/refinement-arc/metawrap*bins/*fa  >> "$1"/prokaryotes/binning/unique_bins/md5_sum ;
cat "$1"/prokaryotes/binning/unique_bins/md5_sum | cut -f1 -d' ' | sort | uniq > "$1"/prokaryotes/binning/unique_bins/md5_unique; 

while read l; 
do 
bininit="$(grep "$l" "$1"/prokaryotes/binning/unique_bins/md5_sum  | head -n1 | cut -f3 -d' ')";
binafter="$1"/prokaryotes/binning/unique_bins/"$lib"-bin."$bin_count".fa; 
cp $bininit $binafter ;
bin_count=$[$bin_count+1];
done < "$1"/prokaryotes/binning/unique_bins/md5_unique ; 
