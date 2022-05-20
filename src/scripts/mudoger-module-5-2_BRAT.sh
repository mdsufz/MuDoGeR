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
output_folder=$project_folder/mapping_results/genome_otu_mapping

# load brat env
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/brat_env

mkdir -p $output_folder

#Copy OTUs based on the results from gOTUpick (5-1) 
mkdir -p $output_folder/otus_fasta

tail -n +2 $project_folder/mapping_results/gOTUpick_results/final_output/bestbins.txt > $output_folder/otus_fasta/aux
while read l;
	do
	otu=$(echo $l | cut -f1 -d " ")
	yes | cp -fr $project_folder/mapping_results/all_bins/$otu $output_folder/otus_fasta
done < $output_folder/otus_fasta/aux
rm -f $output_folder/otus_fasta/aux

if [ "$coverage" = "true" ]; then

	mkdir -p $output_folder/genome_size
	#Calculate genome size
	for bin_path in $output_folder/otus_fasta/*.fa;
	  do  bin="$(echo "$bin_path" | rev | cut -f1 -d'/' | rev )";
	  echo -e "\nCalculating genome size from $bin"
	  cat "$bin_path" | grep -v ">" > $output_folder/genome_size/aux;
	  size="$(wc -c $output_folder/genome_size/aux | cut -f1 -d' ' )";
	  echo ${size} ${bin/.fa/} > "$output_folder/genome_size/${bin/.fa/}.genome-size"
	  rm -f $output_folder/genome_size/aux;

	done
	cat $output_folder/genome_size/*.genome-size > $output_folder/genomes_sizes
fi

#Indexing otus
cd $output_folder/otus_fasta

for bin in *.fa ; 
  do
  echo -e "\nIndexing $bin"
  if [ -f  $output_folder/otus_fasta/${bin/.fa/}.rev.1.bt2 ];
  then echo "-> Bowtie-build already done. Please check here: $output_folder/otus_fasta/${bin/.fa/}.rev.1.bt2"
  else
  echo "-> Running bowtie-build"
  bowtie2-build "$bin" "${bin/.fa/}";
  fi
  
done

cd -

#Prepare reads
mkdir -p $output_folder/merged_reads

cd $output_folder/merged_reads

aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
rm -f $output_folder/merged_reads/total_reads_per_lib.tsv
rm -f $output_folder/merged_reads/avg_reads_len.tsv
for i in $aux; 
	do 
	r1=$project_folder/$i/qc/final_pure_reads_1.fastq; 
	r2=$project_folder/$i/qc/final_pure_reads_2.fastq;
	
	if [ -f  "$output_folder"/merged_reads/$i.fasta ]; then
		echo "-> Reads from $i ready. Please check here: $output_folder/merged_reads/$i.fasta"
  	else
  		echo "-> Merging reads from $i"
  		#Merge reads in pandaseq
  		pandaseq -d rbfkms -T $cores -f "$r1" -r "$r2" -w "$output_folder/merged_reads/$i.fasta"
  	fi
	
	
	#Count number of reads in merged file
	if [ "$relative" = "true" ]; then
		
		#if [ -f  $output_folder/merged_reads/total_reads_per_lib.tsv ]; then
		#	echo "-> Total number of reads from $i counted. Please check here: $output_folder/merged_reads/total_reads_per_lib.tsv"
  		#else
			echo "-> Counting reads from $i"
  			num_reads=`wc -l "$output_folder/merged_reads/$i.fasta"`
  			echo -e "$i\t$num_reads" >> $output_folder/merged_reads/total_reads_per_lib.tsv;
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
			echo -e "$i\t$frag_avg_size" >> $output_folder/merged_reads/avg_reads_len.tsv;
			rm -f aux
		#fi
	
	fi
  
done

cd -

#Run complete if selected
if [ "$complete" = "true" ]; then
	#Create job array for Complete
	cd $output_folder/otus_fasta
	mkdir -p "$output_folder"/map_results_complete/
	
	for d in *.fa; 
		do 
		bin=$output_folder/otus_fasta/${d/.fa/};
		for l in "$output_folder"/merged_reads/*.fasta;
		do 
		output="$output_folder"/map_results_complete/"${d/.fa/}"-LIB-"$(echo $l | rev | cut -f1 -d'/' | rev | sed "s/.fasta/.txt/g")";
		echo "$bin" "$l" "$output" ;
		done;
	done > "$output_folder"/map_list_complete
	
	#Map according to map_list
	while read l;
		do
		bin=$(echo $l | cut -f1 -d " ");
		lib=$(echo $l | cut -f2 -d " ");
		out=$(echo $l | cut -f3 -d " ");
		
		if [[ -f "$out" ]]; then
			
			continue
		else
		
			bowtie2 --threads $cores -x "$bin" -f "$lib" | samtools view --threads $cores  -F 4 | cut -f1 | sort -u | wc -l > "$out"
		fi
		
	done < "$output_folder"/map_list_complete
	cd -
	
	# Create table
	cd "$output_folder"/map_results_complete/
	for d in *.txt;
		do 
		echo -e "$d\t\c" | sed "s/-LIB-/\t/g" | sed "s/.txt//g";
		cat $d;
	done > "$output_folder"/map_results_complete/map_complete_absolute_n_hits_list.tsv
	cd -
	cat "$output_folder"/map_results_complete/map_complete_absolute_n_hits_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$output_folder"/map_results_complete/map_complete_absolute_n_hits_table.tsv
	
	#Calculate coverage and relative abd tables
	if [ "$coverage" = "true" ]; then
		rm -fr "$output_folder"/map_results_complete/map_complete_coverage_list.tsv
		while read l;
		do
			num_hits="$(echo $l | cut -f3 -d " ")";
			lib="$(echo $l | cut -f2 -d " ")";
			bin="$(echo $l | cut -f1 -d " ")";
			frag_size="$(grep -w $lib $output_folder/merged_reads/avg_reads_len.tsv | cut -f2)";
			gen_size="$(grep -w $bin $output_folder/genomes_sizes | cut -f1 -d ' ')";
			hits_times_frag="$(($num_hits*$frag_size))";
			coverage="$(($hits_times_frag/$gen_size))";
			echo "$bin" "$lib" "$coverage" >> "$output_folder"/map_results_complete/map_complete_coverage_list.tsv
		done < "$output_folder"/map_results_complete/map_complete_absolute_n_hits_list.tsv
		
		# transform in cross table 
		cat "$output_folder"/map_results_complete/map_complete_coverage_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$output_folder"/map_results_complete/map_complete_coverage_table.tsv
	fi	
	if [ "$relative" = "true" ]; then
		rm -fr "$output_folder"/map_results_complete/map_complete_relative_abundance_list.tsv
		while read l;
		do 
			num_hits="$(echo $l | cut -f3 -d" ")" ;
			lib="$(echo $l | cut -f2 -d" ")";
			bin="$(echo $l | cut -f1 -d" ")";
			total_n_reads="$(grep -w $lib $output_folder/merged_reads/total_reads_per_lib.tsv | cut -f2 | cut -f1 -d " ")";
			r_abundance="$(  bc -l <<< $num_hits/$total_n_reads)";
			echo "$bin" "$lib" "$r_abundance" >> "$output_folder"/map_results_complete/map_complete_relative_abundance_list.tsv
		done < "$output_folder"/map_results_complete/map_complete_absolute_n_hits_list.tsv 
		# transform in cross table
		cat "$output_folder"/map_results_complete/map_complete_relative_abundance_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$output_folder"/map_results_complete/map_complete_relative_abundance_table.tsv
	fi

fi

#Run reduced if selected

if [ "$reduced" = "true" ]; then
	#Create job array for Reduced
	cd $output_folder/otus_fasta
	mkdir -p "$output_folder"/map_results_reduced/
	
	cat $project_folder/mapping_results/gOTUpick_results/final_output/final_groups_output.csv | grep "*" > "$output_folder"/aux_rep
	rm -f "$output_folder"/map_list_reduced
	touch "$output_folder"/map_list_reduced
	while read l; 
		do
		group=$(echo $l | cut -f2 -d ",");
		rep_bin=$(echo $l | cut -f1 -d ",");
		echo $group;
		if [ $group = "unique" ]; then
			lib=$(echo $rep_bin | cut -f1 -d "-");
			w_bin="$output_folder/otus_fasta/$(echo $rep_bin | sed "s/.fa//g")";
			w_lib="$output_folder/merged_reads/$lib.fasta"
			w_out="$output_folder/map_results_reduced/$(echo $rep_bin | sed "s/.fa//g")-LIB-$lib.txt"
			echo "$w_bin" "$w_lib" "$w_out" >> "$output_folder"/map_list_reduced
		else	
			cat $project_folder/mapping_results/gOTUpick_results/final_output/final_groups_output.csv | grep "$group" > "$output_folder"/aux_group
			while read b;
				do
				lib=$(echo $b | cut -f1 -d "," | cut -f1 -d "-");
				w_bin="$output_folder/otus_fasta/$(echo $rep_bin | sed "s/.fa//g")";
				w_lib="$output_folder/merged_reads/$lib.fasta"
				w_out="$output_folder/map_results_reduced/$(echo $rep_bin | sed "s/.fa//g")-LIB-$lib.txt"
				echo "$w_bin" "$w_lib" "$w_out" >> "$output_folder"/map_list_reduced
			done < "$output_folder"/aux_group
		fi
	done < "$output_folder"/aux_rep
	rm -f "$output_folder"/aux_rep
	rm -f "$output_folder"/aux_group
	
	#Map according to map_list
	while read l;
		do
		bin=$(echo $l | cut -f1 -d " ");
		lib=$(echo $l | cut -f2 -d " ");
		out=$(echo $l | cut -f3 -d " ");
		
		if [[ -f "$out" ]]; then
			
			continue
		else
		
			bowtie2 --threads $cores -x "$bin" -f "$lib" | samtools view --threads $cores  -F 4 | cut -f1 | sort -u | wc -l > "$out"
		fi
		
	done < "$output_folder"/map_list_reduced
	cd -
	
	# Create table
	cd "$output_folder"/map_results_reduced/
	for d in *.txt;
		do 
		echo -e "$d\t\c" | sed "s/-LIB-/\t/g" | sed "s/.txt//g";
		cat $d;
	done > "$output_folder"/map_results_reduced/map_reduced_absolute_n_hits_list.tsv
	cd -
	cat "$output_folder"/map_results_reduced/map_reduced_absolute_n_hits_list.tsv | datamash -sW crosstab 1,2 unique 3 > map_reduced_absolute_n_hits_table.tsv
	
	#Calculate coverage and relative abd tables
	if [ "$coverage" = "true" ]; then
		rm -fr "$output_folder"/map_results_reduced/map_reduced_coverage_list.tsv
		while read l;
		do
			num_hits="$(echo $l | cut -f3 -d" ")";
			lib="$(echo $l | cut -f2 -d" ")";
			bin="$(echo $l | cut -f1 -d" ")";
			frag_size="$(grep -w $lib $output_folder/merged_reads/avg_reads_len.tsv | cut -f2)";
			gen_size="$(grep -w $bin $output_folder/genomes_sizes | cut -f1 -d ' ')";
			hits_times_frag="$(($num_hits*$frag_size))";
			coverage="$(($hits_times_frag/$gen_size))";
			echo "$bin" "$lib" "$coverage" >> "$output_folder"/map_results_reduced/map_reduced_coverage_list.tsv;
		done < "$output_folder"/map_results_reduced/map_reduced_absolute_n_hits_list.tsv 
		
		# transform in cross table 
		cat "$output_folder"/map_results_reduced/map_reduced_coverage_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$output_folder"/map_results_reduced/map_reduced_coverage_table.tsv
	fi
	if [ "$relative" = "true" ]; then
		rm -fr "$output_folder"/map_results_reduced/map_reduced_relative_abundance_list.tsv
		while read l;
		do 
			num_hits="$(echo $l | cut -f3 -d" ")" ;
			lib="$(echo $l | cut -f2 -d" ")";
			bin="$(echo $l | cut -f1 -d" ")";
			total_n_reads="$(grep -w $lib $output_folder/merged_reads/total_reads_per_lib.tsv | cut -f2 | cut -f1 -d " ")";
			r_abundance="$(  bc -l <<< $num_hits/$total_n_reads)";
			echo "$bin" "$lib" "$r_abundance" >> "$output_folder"/map_results_reduced/map_reduced_relative_abundance_list.tsv
		done < "$output_folder"/map_results_reduced/map_reduced_absolute_n_hits_list.tsv 
		# transform in cross table
		cat "$output_folder"/map_results_reduced/map_reduced_relative_abundance_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$output_folder"/map_results_reduced/map_reduced_relative_abundance_table.tsv
	fi

fi

#deactivate env
conda deactivate




