#!/bin/bash

# This bash script calculates the BRAT

# loading conda environment
echo '------- START MODULE 5-2-2 UViGs BRAT'

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
uvigs_output_folder=$project_folder/mapping_results/uvigs_mapping
merged_reads_folder=$project_folder/mapping_results/merged_reads

# load brat env
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/brat_env

mkdir -p $uvigs_output_folder

#### Copy seq files ####
mkdir -p $uvigs_output_folder/uvigs_fasta

#Copy recovered UViGs 

aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
for i in $aux; 
	do
	for f in $project_folder/$i/viruses/final_outputs/only_uvigs_seq/*.fa
		do
		file=` echo $f | rev | cut -f1 -d'/' | rev`
		yes | cp $f $uvigs_output_folder/uvigs_fasta/"$i"-"$file"
	done
done

#### Calculate MAGs/UViGs/eMABs sizes ###

if [ "$coverage" = "true" ]; then

	## UViGs sizes ##
	mkdir -p $uvigs_output_folder/uvigs_sizes
	
	for uvigs_path in $uvigs_output_folder/uvigs_fasta/*.fa;
	  do  uvig="$(echo "$uvigs_path" | rev | cut -f1 -d'/' | rev )";
	  echo -e "\nCalculating genome size from $uvig"
	  cat "$uvigs_path" | grep -v ">" > $uvigs_output_folder/uvigs_sizes/aux;
	  size="$(wc -c $uvigs_output_folder/uvigs_sizes/aux | cut -f1 -d' ' )";
	  echo ${size} ${uvig/.fa/} > "$uvigs_output_folder/uvigs_sizes/${uvig/.fa/}.uvig-size"
	  rm -f $uvigs_output_folder/uvigs_sizes/aux;

	done
	cat $uvigs_output_folder/uvigs_sizes/*.uvig-size > $uvigs_output_folder/uvigs_sizes.txt
		
fi

#### Indexing OTUs/UViGs/eMABs ####

## Indexing UViGs
cd $uvigs_output_folder/uvigs_fasta/

for uvig in *.fa ; 
  do
  echo -e "\nIndexing $uvig"
  if [ -f  $uvigs_output_folder/uvigs_fasta/${uvig/.fa/}.rev.1.bt2 ];
  then echo "-> Bowtie-build already done. Please check here: $uvigs_output_folder/uvigs_fasta/${uvig/.fa/}.rev.1.bt2"
  else
  echo "-> Running bowtie-build"
  bowtie2-build "$uvig" "${uvig/.fa/}";
  fi
  
done

cd -

#Run complete if selected
if [ "$complete" = "true" ]; then
	#Create job array for Complete
	cd $uvigs_output_folder/uvigs_fasta
	mkdir -p "$uvigs_output_folder"/map_results_complete/
	
	for d in *.fa; 
		do 
		bin=$uvigs_output_folder/uvigs_fasta/${d/.fa/};
		for l in $merged_reads_folder/*.fasta;
		do 
		output="$uvigs_output_folder"/map_results_complete/"${d/.fa/}"-LIB-"$(echo $l | rev | cut -f1 -d'/' | rev | sed "s/.fasta/.txt/g")";
		echo "$bin" "$l" "$output" ;
		done;
	done > "$uvigs_output_folder"/map_list_complete
	
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
		
	done < "$uvigs_output_folder"/map_list_complete
	cd -
	rm -f "$uvigs_output_folder"/map_list_complete
  
	# Create table
	cd "$uvigs_output_folder"/map_results_complete/
	for d in *.txt;
		do 
		echo -e "$d\t\c" | sed "s/-LIB-/\t/g" | sed "s/.txt//g";
		cat $d;
	done > "$uvigs_output_folder"/map_results_complete/map_complete_absolute_n_hits_list.tsv
	cd -
	cat "$uvigs_output_folder"/map_results_complete/map_complete_absolute_n_hits_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$uvigs_output_folder"/map_results_complete/map_complete_absolute_n_hits_table.tsv
	
	#Calculate coverage and relative abd tables
	if [ "$coverage" = "true" ]; then
		rm -fr "$uvigs_output_folder"/map_results_complete/map_complete_coverage_list.tsv
		while read l;
		do
			num_hits="$(echo $l | cut -f3 -d " ")";
			lib="$(echo $l | cut -f2 -d " ")";
			bin="$(echo $l | cut -f1 -d " ")";
			frag_size="$(grep -w $lib $merged_reads_folder/avg_reads_len.tsv | cut -f2)";
			gen_size="$(grep -w $bin $uvigs_output_folder/uvigs_sizes.txt | cut -f1 -d ' ')";
			hits_times_frag="$(($num_hits*$frag_size))";
			coverage="$(($hits_times_frag/$gen_size))";
			echo "$bin" "$lib" "$coverage" >> "$uvigs_output_folder"/map_results_complete/map_complete_coverage_list.tsv
		done < "$uvigs_output_folder"/map_results_complete/map_complete_absolute_n_hits_list.tsv
		
		# transform in cross table 
		cat "$uvigs_output_folder"/map_results_complete/map_complete_coverage_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$uvigs_output_folder"/map_results_complete/map_complete_coverage_table.tsv
	fi
  
	if [ "$relative" = "true" ]; then
		rm -fr "$uvigs_output_folder"/map_results_complete/map_complete_relative_abundance_list.tsv
		while read l;
		do 
			num_hits="$(echo $l | cut -f3 -d" ")" ;
			lib="$(echo $l | cut -f2 -d" ")";
			bin="$(echo $l | cut -f1 -d" ")";
			total_n_reads="$(grep -w $lib $merged_reads_folder/total_reads_per_lib.tsv | cut -f2 | cut -f1 -d " ")";
			r_abundance="$(  bc -l <<< $num_hits/$total_n_reads)";
			echo "$bin" "$lib" "$r_abundance" >> "$uvigs_output_folder"/map_results_complete/map_complete_relative_abundance_list.tsv
		done < "$uvigs_output_folder"/map_results_complete/map_complete_absolute_n_hits_list.tsv 
		# transform in cross table
		cat "$uvigs_output_folder"/map_results_complete/map_complete_relative_abundance_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$uvigs_output_folder"/map_results_complete/map_complete_relative_abundance_table.tsv
	fi
mkdir -p $uvigs_output_folder/map_final_tables_complete
mv -f "$uvigs_output_folder"/map_results_complete/map_complete_* $uvigs_output_folder/map_final_tables_complete
fi

### Run reduced if selected ###

if [ "$reduced" = "true" ]; then
	#Create job array for Reduced
	cd $uvigs_output_folder/uvigs_fasta
	mkdir -p "$uvigs_output_folder"/map_results_reduced/
	
	rm -f "$uvigs_output_folder"/map_list_reduced
	touch "$uvigs_output_folder"/map_list_reduced
	
	for f in $uvigs_output_folder/uvigs_fasta/*.fa;
		do 
		lib=`echo $f | rev | cut -f1 -d'/' | rev | cut -f1 -d'-'`;
		bin=`echo $f | rev | cut -f1 -d'/' | rev | sed "s/.fa//g"`;
		w_bin="$uvigs_output_folder/uvigs_fasta/$bin";
		w_lib="$merged_reads_folder/$lib.fasta"
		w_out="$uvigs_output_folder/map_results_reduced/$bin-LIB-$lib.txt"
		echo "$w_bin" "$w_lib" "$w_out" >> "$uvigs_output_folder"/map_list_reduced
	done
	
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
		
	done < "$uvigs_output_folder"/map_list_reduced
	cd -
	rm -fr "$uvigs_output_folder"/map_list_reduced
	# Create table
	cd "$uvigs_output_folder"/map_results_reduced/
	for d in *.txt;
		do 
		echo -e "$d\t\c" | sed "s/-LIB-/\t/g" | sed "s/.txt//g";
		cat $d;
	done > "$uvigs_output_folder"/map_results_reduced/map_reduced_absolute_n_hits_list.tsv
	cd -
	cat "$uvigs_output_folder"/map_results_reduced/map_reduced_absolute_n_hits_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$uvigs_output_folder"/map_results_reduced/map_reduced_absolute_n_hits_table.tsv
	
	#Calculate coverage and relative abd tables
	if [ "$coverage" = "true" ]; then
		rm -fr "$uvigs_output_folder"/map_results_reduced/map_reduced_coverage_list.tsv
		while read l;
		do
			num_hits="$(echo $l | cut -f3 -d" ")";
			lib="$(echo $l | cut -f2 -d" ")";
			bin="$(echo $l | cut -f1 -d" ")";
			frag_size="$(grep -w $lib $merged_reads_folder/avg_reads_len.tsv | cut -f2)";
			gen_size="$(grep -w $bin $uvigs_output_folder/uvigs_sizes.txt | cut -f1 -d ' ')";
			hits_times_frag="$(($num_hits*$frag_size))";
			coverage="$(($hits_times_frag/$gen_size))";
			echo "$bin" "$lib" "$coverage" >> "$uvigs_output_folder"/map_results_reduced/map_reduced_coverage_list.tsv;
		done < "$uvigs_output_folder"/map_results_reduced/map_reduced_absolute_n_hits_list.tsv 
		
		# transform in cross table 
		cat "$uvigs_output_folder"/map_results_reduced/map_reduced_coverage_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$uvigs_output_folder"/map_results_reduced/map_reduced_coverage_table.tsv
	fi
	if [ "$relative" = "true" ]; then
		rm -fr "$uvigs_output_folder"/map_results_reduced/map_reduced_relative_abundance_list.tsv
		while read l;
		do 
			num_hits="$(echo $l | cut -f3 -d" ")" ;
			lib="$(echo $l | cut -f2 -d" ")";
			bin="$(echo $l | cut -f1 -d" ")";
			total_n_reads="$(grep -w $lib $merged_reads_folder/total_reads_per_lib.tsv | cut -f2 | cut -f1 -d " ")";
			r_abundance="$(  bc -l <<< $num_hits/$total_n_reads)";
			echo "$bin" "$lib" "$r_abundance" >> "$uvigs_output_folder"/map_results_reduced/map_reduced_relative_abundance_list.tsv
		done < "$uvigs_output_folder"/map_results_reduced/map_reduced_absolute_n_hits_list.tsv 
		# transform in cross table
		cat "$uvigs_output_folder"/map_results_reduced/map_reduced_relative_abundance_list.tsv | datamash -sW crosstab 1,2 unique 3 > "$uvigs_output_folder"/map_results_reduced/map_reduced_relative_abundance_table.tsv
	fi
mkdir -p $uvigs_output_folder/map_final_tables_reduced
mv "$uvigs_output_folder"/map_results_reduced/map_reduced_* $uvigs_output_folder/map_final_tables_reduced
fi

#deactivate env
conda deactivate

