#!/bin/bash

### VIRAL INVESTIGATION - VIRSORTER, VIRFINDER, VIBRANT And STAMPEDE DEREPLICATION
echo -e '\n\n------- START MODULE 3-1 VIRAL INVESTIGATION \n\n'

# arguments declaration
log="log_vir"                   # definition of path to log file
assembly=$1                     # assembly file
viruses_folder=$2               # master virus folder
output_folder=$2/investigation  # output folder to be created inside master output folder
num_cores=$3                    # number of threads

conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

lib=`echo $output_folder | rev | cut -f4 -d'/' | rev `

############ VIBRANT #################

if [ -f "$output_folder"/vibrant/VIBRANT_final_assembly/VIBRANT_phages_final_assembly/final_assembly.phages_combined.fna ];
then echo -e '\n----> Vibrant investigation is done\n'
else
echo "-----> START VIBRANT (1/4)"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/vibrant_env
conda_vib="$(echo $PATH | cut -f1 -d':')"
VIBRANT_FILES="$(echo $PATH | sed "s/:/\n/g" | grep vibrant_env  | sed "s/bin/share\/vibrant-1.2.0\/files/g")"
VIBRANT_run.py -i $assembly -folder "$output_folder"/vibrant -t $num_cores -m "$VIBRANT_FILES" -d $DATABASES_LOCATION/vibrant
# fetch results
cat "$output_folder"/vibrant/VIBRANT_final_assembly/VIBRANT_phages_final_assembly/final_assembly.phages_combined.fna | 
grep ">" | sed "s/_fragment_1//g;s/>//g"   > "$output_folder"/vibrant_filtered_data.txt
conda deactivate
echo "-----> END VIBRANT (1/4)"
fi

######### VIRFINDER        
dependencies="$(echo $PATH | cut -f1 -d':' | sed "s/bin/dependencies/g")"
if [ -f "$output_folder"/virfinder/virfinder_output.tsv ];
then echo -e '\n----> Virfinder investigation is done\n'
else
echo "-----> START VIRFINDER (2/4)"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/virfinder_env
assembly_whole_path="$(realpath "$assembly")"
output_file="$(realpath "$output_folder"/virfinder)"/virfinder_output.tsv
mkdir -p "$output_folder"/virfinder
Rscript $dependencies/vir_virfinder_script.r "$output_folder"/virfinder "$assembly_whole_path" "$output_file"
# fetch results
cat "$output_folder"/virfinder/virfinder_output.tsv | awk -F'\t' '{ if ( $4 <= 0.01) print }' | 
awk -F'_' '{ if ( $4 >= 1000) print  }' | cut -f2 | sed "s/\"//g" > "$output_folder"/virfinder_filtered_data.txt
conda deactivate
echo "-----> END VIRFINDER (2/4)"
fi


######### VIRSORTER
if [ -f "$output_folder"/virsorter/final-viral-combined.fa ] ;
then echo -e '\n----> Virsoter investigation is done\n'
else
echo "-----> START VIRSORTER (3/4)"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/virsorter2_env
virsorter run all -i "$assembly" -w "$output_folder"/virsorter -j "$num_cores"
# fetch results
cat "$output_folder"/virsorter/final-viral-combined.fa  | grep ">" | sed "s/_fragment_1//g;s/>//g" | 
cut -f1 -d "|" > "$output_folder"/virsorter2_filtered_data.txt
conda deactivate
echo "-----> END VIRSORTER (3/4)"
fi


####### GET RESULTS AND PUT THEM TOGETHER
mkdir -p "$output_folder"/dereplication
cat $output_folder/*txt | sort -u >  "$output_folder"/dereplication/viral_unique_contigs


### PICK FASTA SEQUENCES OF VIRAL CONTIGS FOR FURTHER DEREPLICATION

if [ -f "$output_folder"/dereplication/uvigs.fa ];
then echo -e '\n----> Viral dereplication is done\n'
else
echo "-----> START DEREPLICATION (4/4)"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/extract_env
python "$MUDOGER_DEPENDENCIES_PATH"/extract_fa.py "$output_folder"/dereplication/viral_unique_contigs "$assembly" "$output_folder"/dereplication/uvigs.fa
conda deactivate

fi

##### RUN DEREPLICATION
if [ ! -f "$output_folder"/dereplication/uvigs_95-70.fna ];
then
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/stampede_clustergenomes_env
Cluster_genomes.pl -f "$output_folder"/dereplication/uvigs.fa  -c 70 -i 95
mkdir -p "$viruses_folder"/final_outputs/putative_viral_contigs
python3 "$MUDOGER_DEPENDENCIES_PATH"/split-all-seq.py "$output_folder"/dereplication/uvigs_95-70.fna "$viruses_folder"/final_outputs/putative_viral_contigs/"$lib"_putative_viral_contig
conda deactivate
echo "-----> END DEREPLICATION (4/4)"
else
echo "-> Dereplication is done";
fi

