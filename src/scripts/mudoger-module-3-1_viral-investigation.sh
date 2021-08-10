#!/bin/bash

### VIRAL INVESTIGATION - VIRSORTER, VIRFINDER, VIBRANT And STAMPEDE DEREPLICATION

# arguments declaration
log="log_vir"                      # definition of path to QC log       
assembly=$1
#forward_library=$2              # forward library path  # commented because reads are not used
#reverse_library=$3              # reverse library path   # and thus not necessary
output_folder=$2                # output folder to be created inside master output folder
num_cores=$3                     # number of threads


############ VIBRANT
conda activate vibrant-env
python VIBRANT_run.py -i $assembly -folder "$output_folder"/vibrant -t $num_cores
# fetch results
cat "$output_folder"/vibrant/VIBRANT_final_assembly/VIBRANT_phages_final_assembly/final_assembly.phages_combined.fna | 
grep ">" | sed "s/_fragment_1//g;s/>//g"   > "$output_folder"/VIBRANT_filtered_data.txt
conda deactivate

######### VIRFINDER
conda activate virfinder-env
mkdir -p "$output_folder"/virfinder
Rscript MuDoGeR/tools/virfinder_script.r "$output_folder"/virfinder "$assembly" "$output_folder"/virfinder/virfinder_output.tsv
# fetch results
cat "$output_folder"/virfinder/virfinder_output.tsv | awk -F'\t' '{ if ( $4 <= 0.01) print }' | 
awk -F'_' '{ if ( $4 >= 1000) print  }' | cut -f2 | sed "s/\"//g" > "$output_folder"/virfinder/virfinder_filtered_data.txt
conda deactivate

######### VIRSORTER
conda activate virsorter2-env
virsorter run all -i "$assembly" -w "$output_folder"/virsorter -j "$num_cores"
# fetch results
cat "$output_folder"/virsorter/final-viral-combined.fa  | grep ">" | sed "s/_fragment_1//g;s/>//g" | 
cut -f1 -d "|" > "$output_folder"/virsorter/virsorter2_filtered_data.txt
conda deactivate

####### GET RESULTS AND PUT THEM TOGETHER
mkdir -p "$output_folder"/dereplication
cat "$output_folder"/VIBRANT_filtered_data.txt \
"$output_folder"/virfinder/virfinder_filtered_data.txt \ 
"$output_folder"/virsorter/virsorter2_filtered_data.txt | sort -u >  "$output_folder"/dereplication/viral_unique_contigs

### PICK FASTA SEQUENCES OF VIRAL CONTIGS FOR FURTHER DEREPLICATION
conda activate extract-env
python MuDoGeR/tools/extract_fa.py "$output_folder"/dereplication/viral_unique_contigs "$assembly" "$output_folder"/dereplication/uvigs.fa
conda deactivate

##### RUN DEREPLICATION
conda activate stampede-clustergenomes-env
Cluster_genomes.pl -f "$output_folder"/dereplication/uvigs.fa  -c 70 -i 95
conda deactivate
 
