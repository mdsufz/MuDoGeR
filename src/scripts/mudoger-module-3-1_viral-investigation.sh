
### VIRAL INVESTIGATION - VIRSORTER, VIRFINDER, VIBRANT And STAMPEDE DEREPLICATION

# arguments declaration
log="log_vir"                      # definition of path to QC log       
assembly=$1
forward_library = $2              # forward library path
reverse_library = $3              # reverse library path
output_folder = $4                # output folder to be created inside master output folder
num_cores = $5                     # number of threads


############ VIBRANT
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 Python/3.6.6
source /data/msb/tools/vibrant/env_vibrant_v1.2.1/bin/activate
python /data/msb/tools/vibrant/VIBRANT/VIBRANT_run.py -i $assembly -folder "$output_folder"/vibrant -t $num_cores
# fetch results
cat "$output_folder"/vibrant/VIBRANT_final_assembly/VIBRANT_phages_final_assembly/final_assembly.phages_combined.fna | 
grep ">" | sed "s/_fragment_1//g;s/>//g"   > "$output_folder"/VIBRANT_filtered_data.txt


######### VIRFINDER 
mkdir -p "$output_folder"/virfinder
Rscript /data/msb/tools/virfinder/virfinder_script.r "$output_folder"/virfinder "$assembly" "$output_folder"/virfinder/virfinder_output.tsv
# fetch results
cat "$output_folder"/virfinder/virfinder_output.tsv | awk -F'\t' '{ if ( $4 <= 0.01) print }' | 
awk -F'_' '{ if ( $4 >= 1000) print  }' | cut -f2 | sed "s/\"//g" > "$output_folder"/virfinder/virfinder_filtered_data.txt

######### VIRSORTER
source activate /data/msb/tools/virsorter2/virsorter2-conda
virsorter run all -i "$assembly" -w "$output_folder"/virsorter -j "$num_cores"
# fetch results
cat "$output_folder"/virsorter2/final-viral-combined.fa  | grep ">" | sed "s/_fragment_1//g;s/>//g" | 
cut -f1 -d "|" > "$output_folder"/virsorter/virsorter2_filtered_data.txt

####### GET RESULTS AND PUT THEM TOGETHER
mkdir -p "$output_folder"/dereplication
cat "$output_folder"/VIBRANT_filtered_data.txt \
"$output_folder"/virfinder/virfinder_filtered_data.txt \ 
"$output_folder"/virsorter/virsorter2_filtered_data.txt | sort -u >  "$output_folder"/dereplication/viral_unique_contigs

### EXTRACT 
python extract_fa.py "$output_folder"/dereplication/viral_unique_contigs "$assembly" "$output_folder"/dereplication/uvigs.fa

##### RUN DEREPLICATION

 
