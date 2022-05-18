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
meta_table=$2
cores=$3
absolute=$4;
coverage=$5;
relative=$6;
reduced=$7
complete=$8

#Process inputs
output_folder=$project_folder/mapping_results/genome_bins_mapping

# load brat env
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/brat_env

mkdir -p $output_folder

#Copy bins and Calculate genome size
mkdir -p $output_folder/genome_size

for bin_path in $project_folder/mapping_results/all_bins/*.fa;
  do  bin="$(echo "$bin_path" | rev | cut -f1 -d'/' | rev )";
  echo -e "\nCalculating genome size from $bin"
  cat "$bin_path" | grep -v ">" > $output_folder/genome_size/aux;
  size="$(wc -c $output_folder/genome_size/aux | cut -f1 -d' ' )";
  echo ${size} ${bin/.fa/} > "$output_folder/genome_size/${bin/.fa/}.genome-size"
  rm -f $output_folder/genome_size/aux;
  
done
cat $output_folder/genome_size/*.genome-size > $output_folder/genomes_sizes

#Indexing bins
cd $project_folder/mapping_results/all_bins/

for bin in *.fa ; 
  do
  echo -e "\nIndexing $bin"
  if [ -f  $project_folder/mapping_results/all_bins/${bin/.fa/}.rev.1.bt2 ];
  then echo "-> Bowtie-build already done. Please check here: $project_folder/mapping_results/all_bins/${bin/.fa/}.rev.1.bt2"
  else
  echo "-> Running bowtie-build"
  bowtie2-build "$bin" "${bin/.fa/}";
  fi
  
done

cd -

#Merge Pair-end reads
mkdir -p $output_folder/merged_reads

cd $output_folder/merged_reads

aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
for i in $aux; 
	do 
	r1=$project_folder/$i/qc/final_pure_reads_1.fastq; 
	r2=$project_folder/$i/qc/final_pure_reads_2.fastq;
	
	#Merge reads in pandaseq
  	pandaseq -f "$r1" -r "$r2" -w "$output_folder/merged_reads/$i.fasta"
	
	#Count number of reads in merged file
  	echo "Calculating number of reads from $i"
  	num_reads=`wc -l "$output_folder/merged_reads/$i.fasta"`
  	echo -e "$i\t$num_reads" >> $output_folder/merged_reads/total_reads_per_lib.tsv;
	
	#Calculate average read size in lib
	cat $i.fasta | grep -v ">" > aux;
	seqs_num="$(wc -l aux | cut -f1 -d' ' )" ;
	size="$(wc -c aux | cut -f1 -d' ' )";
	frag_avg_size="$((size/seqs_num))";
	echo -e "$i\t$frag_avg_size" >> $output_folder/merged_reads/avg_reads_len.tsv;
	rm -f aux
  
done

cd -

#Create job array for Complete
cd $project_folder/mapping_results/all_bins/

for d in *.fa; 
	do 
	bin=$project_folder/mapping_results/all_bins/${d/.fa/};
	for l in "$output_folder"/merged_reads/*.fasta;
	do 
	output="$output_folder"/mappings/"${d/.fa/}"-LIB-"$(echo $l | rev | cut -f1 -d'/' | rev | sed "s/.fasta/.txt/g")";
	echo "$bin" "$l" "$output" ;
	done;
done > "$output_folder"/map_list

#deactivate env
conda deactivate




