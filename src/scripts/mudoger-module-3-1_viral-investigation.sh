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
echo "-----> STARTING VIBRANT (1/4)"
if [ -f "$output_folder"/vibrant/VIBRANT_final_assembly/VIBRANT_phages_final_assembly/final_assembly.phages_combined.fna ];
then echo '-> Vibrant investigation is done.'
else
conda activate vibrant-env
conda_vib="$(echo $PATH | cut -f1 -d':')"
echo VIBRANT_run.py -i $assembly -folder "$output_folder"/vibrant -t $num_cores -m '$conda_vib'/files -d '$conda_vib'/databases
# fetch results
cat "$output_folder"/vibrant/VIBRANT_final_assembly/VIBRANT_phages_final_assembly/final_assembly.phages_combined.fna | 
grep ">" | sed "s/_fragment_1//g;s/>//g"   > "$output_folder"/vibrant_filtered_data.txt
conda deactivate
echo "-----> END VIBRANT (1/4)"
fi

exit 0

######### VIRFINDER        
echo "-----> STARTING VIRFINDER (2/4)"
conda activate virfinder-env
assembly_whole_path="$(realpath "$assembly")"
output_file="$(realpath "$output_folder"/virfinder)"/virfinder_output.tsv
mkdir -p "$output_folder"/virfinder
Rscript MuDoGeR/tools/vir_virfinder_script.r "$output_folder"/virfinder "$assembly_whole_path" "$output_file"
# fetch results
cat "$output_folder"/virfinder/virfinder_output.tsv | awk -F'\t' '{ if ( $4 <= 0.01) print }' | 
awk -F'_' '{ if ( $4 >= 1000) print  }' | cut -f2 | sed "s/\"//g" > "$output_folder"/virfinder_filtered_data.txt
conda deactivate
echo "-----> END VIRFINDER (2/4)"


######### VIRSORTER
echo "-----> STARTING VIRSORTER (3/4)"
conda activate virsorter2-env
virsorter run all -i "$assembly" -w "$output_folder"/virsorter -j "$num_cores"
# fetch results
cat "$output_folder"/virsorter/final-viral-combined.fa  | grep ">" | sed "s/_fragment_1//g;s/>//g" | 
cut -f1 -d "|" > "$output_folder"/virsorter2_filtered_data.txt
conda deactivate
echo "-----> END VIRSORTER (3/4)"

####### GET RESULTS AND PUT THEM TOGETHER
mkdir -p "$output_folder"/dereplication
cat $output_folder/*txt | sort -u >  "$output_folder"/dereplication/viral_unique_contigs


### PICK FASTA SEQUENCES OF VIRAL CONTIGS FOR FURTHER DEREPLICATION
echo "-----> STARTING DEREPLICATION (4/4)"
conda activate extract-env
python MuDoGeR/tools/extract_fa.py "$output_folder"/dereplication/viral_unique_contigs "$assembly" "$output_folder"/dereplication/uvigs.fa
conda deactivate

##### RUN DEREPLICATION
conda activate stampede-clustergenomes-env
Cluster_genomes.pl -f "$output_folder"/dereplication/uvigs.fa  -c 70 -i 95
conda deactivate
echo "-----> END DEREPLICATION (4/4)"
