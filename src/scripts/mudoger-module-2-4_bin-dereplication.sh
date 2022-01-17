#!/bin/bash

# 4 BIN REDUNANCY REMOVAL
echo '------- START MODULE 2-4 BIN REDUNDANCY REMOVAL'
lib_folder=$1;
bin_count=0;

lib_name=$(basename $lib_folder);

md5sum "$lib_folder"/prokaryotes/binning/refinement-bac/metawrap*bins/*fa  >  "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum ;
md5sum "$lib_folder"/prokaryotes/binning/refinement-arc/metawrap*bins/*fa  >> "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum ;
cat "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum | cut -f1 -d' ' | sort | uniq > "$lib_folder"/prokaryotes/binning/unique_bins/md5_unique; 

while read l; 
do 
bininit="$(grep "$l" "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum  | head -n1 | cut -f3 -d' ')";
binafter="$lib_folder"/prokaryotes/binning/unique_bins/"$lib_name"-bin."$bin_count".fa; 
cp $bininit $binafter ;
bin_count=$[$bin_count+1];
done < "$lib_folder"/prokaryotes/binning/unique_bins/md5_unique ; 

rm -f "$lib_folder"/prokaryotes/binning/unique_bins/md5_unique
rm -f "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum
