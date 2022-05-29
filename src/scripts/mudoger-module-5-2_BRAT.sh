#!/bin/bash

# This bash script calculates the BRAT

# loading conda environment
echo '------- START MODULE 5-2 BRAT'

conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

#Input data
project_folder=$1
metadata_table=$2
cores=$3
absolute=$4;
coverage=$5;
relative=$6;
reduced=$7
complete=$8

#Process inputs
output_folder=$project_folder/mapping_results/pmags_otu_mapping
uvigs_output_folder=$project_folder/mapping_results/uvigs_mapping
emabs_output_folder=$project_folder/mapping_results/euk_mabs_mapping

# load brat env
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/brat_env

mkdir -p $output_folder
mkdir -p $uvigs_output_folder
mkdir -p $emabs_output_folder


#### Prepare reads ####
mkdir -p $project_folder/mapping_results/merged_reads
merged_reads_folder=$project_folder/mapping_results/merged_reads

cd $merged_reads_folder

aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
rm -f $merged_reads_folder/total_reads_per_lib.tsv
rm -f $merged_reads_folder/avg_reads_len.tsv
for i in $aux; 
	do 
	r1=$project_folder/$i/qc/final_pure_reads_1.fastq; 
	r2=$project_folder/$i/qc/final_pure_reads_2.fastq;
	
	if [ -f  $merged_reads_folder/$i.fasta ]; then
		echo "-> Reads from $i ready. Please check here: $merged_reads_folder/$i.fasta"
  	else
  		echo "-> Merging reads from $i"
  		#Merge reads in pandaseq
  		pandaseq -d rbfkms -T $cores -f "$r1" -r "$r2" -w "$merged_reads_folder/$i.fasta"
  	fi
	
	
	#Count number of reads in merged file
	if [ "$relative" = "true" ]; then
		
		#if [ -f  $output_folder/merged_reads/total_reads_per_lib.tsv ]; then
		#	echo "-> Total number of reads from $i counted. Please check here: $output_folder/merged_reads/total_reads_per_lib.tsv"
  		#else
			echo "-> Counting reads from $i"
  			num_reads=`wc -l "$merged_reads_folder/$i.fasta"`
  			echo -e "$i\t$num_reads" >> $merged_reads_folder/total_reads_per_lib.tsv;
		#fi
	fi
	
	#Calculate average read size in lib
	if [ "$coverage" = "true" ]; then
	
		#if [ -f  $output_folder/merged_reads/avg_reads_len.tsv ]; then
		#	echo "-> Average read size from $i calculated. Please check here: $output_folder/merged_reads/avg_reads_len.tsv "
		#else
			echo "-> Calculating average read size from $i"
			cat $i.fasta | grep -v ">" > aux;
			seqs_num="$(wc -l aux | cut -f1 -d' ' )";
			size="$(wc -c aux | cut -f1 -d' ' )";
			frag_avg_size="$((size/seqs_num))";
			echo -e "$i\t$frag_avg_size" >> $merged_reads_folder/avg_reads_len.tsv;
			rm -f aux
		#fi
	
	fi
  
done

cd -

#     Perform BRAT      submodule 5-2-1
if [ -f $project_folder/mapping_results/gOTUpick_results/final_output/bestbins.txt ]; then

	echo -e "\n BRAT pMAGs STARTED"
	bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-5-2-1_BRAT_pmags.sh $project_folder \
										$metadata_table \
										$cores \
										$absolute \
										$coverage \
										$relative \
										$reduced \
										$complete
	echo -e "\n BRAT pMAGs DONE"
	
	
#     Perform BRAT      submodule 5-2-2
if [ -f $project_folder/mapping_results/gOTUpick_results/final_output/bestbins.txt ]; then

	echo -e "\n BRAT UViGs STARTED"
	bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-5-2-2_BRAT_uvigs.sh $project_folder \
										$metadata_table \
										$cores \
										$absolute \
										$coverage \
										$relative \
										$reduced \
										$complete
	echo -e "\n BRAT UViGs DONE"
	
	
#     Perform BRAT      submodule 5-2-3
#if [ -f $project_folder/mapping_results/gOTUpick_results/final_output/bestbins.txt ]; then

	#echo -e "\n BRAT eMABs STARTED"
	#bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-5-2-3_BRAT_emabs.sh $project_folder \
	#									$metadata_table \
	#									$cores \
	#									$absolute \
	#									$coverage \
	#									$relative \
	#									$reduced \
	#									$complete
	#echo -e "\n BRAT eMABs DONE"


#deactivate env
conda deactivate




