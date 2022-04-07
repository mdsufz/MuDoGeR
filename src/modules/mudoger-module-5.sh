#!/bin/bash

#BRAT: a tool to calculate relative abundance of fasta files in fastq files
#MERGE RESULTS and run BRAT
#Running bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/$module_script -o $output_folder -t $num_cores $brat_type;
#The user should have initial data to use BRAT: output from GTDB-Tk, output from CheckM, output from BBMap, the recovered bins/FASTA in fasta format, and Quality Controlled Pair End reads.

help_message() {
		echo ""
		echo "Example usage: BRAT.sh  -f path/to/fasta_files_folder -r path/to/reads_fastq_folder -o path/to/outputdir --reduced --coverage"
		echo ""
		echo "You must select one BRAT type: --reduced  or --complete BRAT"
		echo ""
		echo "You must select at least one output option: --absolute-values and/or --coverage and/or --relative-abundance"
		echo ""
		echo "Available options:"
		echo ""
		echo "-h, --help"
		echo "	display this help page"
		echo ""
		echo "--meta <value>"
		echo "	absolute path to metadata file as explained in MuDoGeR github documentation"
		echo "-f <value>"
		echo "	absolute path to a directory where fasta files are are stored (.fa)"
		echo ""
		echo "-r <value>"
		echo "	absolute path to a directory where reads in fastq format are stored (.fastq)"
		echo ""
		echo "-o <value>"
		echo "	absolute path to a directory where the user wants to store the output"
		echo ""
		echo ""
		echo "-t <value>"
		echo "	number of threads/cores [default: 1]"
		echo ""
		echo "BRAT types:"
		echo ""
		echo "--reduced"	
		echo "	run BRAT of the fasta files only to the origin’s read. For prokaryoten, it will be the representative OTU to the reads of all other OTU inside the same species’ group"
		echo ""
		echo "--complete"	
		echo "	run BRAT of the fasta files to all reads."
		echo ""
		echo "Output types:"
		echo ""
		echo "--absolute-values"	
		echo "	calculate absolute number of reads mapped."
		echo ""
		echo "--coverage"	
		echo "	calculate coverage."
		echo ""
		echo "--relative-abundance"	
		echo "	calculate relative abundance."
		echo "";
				}

#Default parameters
#"WORKDIR" is the path to the output directory.
#"FASTA" is the path to the directory where FASTA/bins are found. 
#"READS" is the path to the directory where samples are found. 

#WORKDIR=None; FASTA=None; READS=None; THREAD=1; reduced=true; complete=false; absolute=false; coverage=false; relative=false

OPTS='getopt -o f:r:o:t --long reduced:,complete:,absolute-values:,coverage:,relative-abundance:,help -- "$@" '

if [ $? -ne 0 ]; then echo "Warning: Something unexpected happened" help_message; exit 1; fi

while true 
do
	case $1 in
		#-f) FASTA="$2"; shift 2;;
		#-r) READS="$2"; shift 2;;
		-o) WORKDIR="$2"; shift 2;;
		-t) THREAD="$2"; shift 2;;
		--meta) metadata_table=$2; shift 2;;
		--reduced) reduced=true; shift 1;; 
		--complete) complete=true; shift 1;; 
		--absolute-values) absolute=true; shift 1;; 
		--coverage) coverage=true; shift 1;; 
		--relative-abundance) relative=true; shift 1;; 
		-h | --help) help_message; exit 1; shift 1;;
		*) break;;
	esac
done

#Checking if output path are entered correctly
if [ "$WORKDIR" = "None" ] || [ "$metadata_table" = "None" ] ; then 
	echo "Non-optional parameter for BRAT type was not entered. Please select output path -o path/to/output/folder and metadata table"
	help_message; exit 1
fi

#Checking if all non-optional parameters are entered correctly
#if [ "$FASTA" = "None" ]  || [ "$READS" = "None" ] ; then 
#	echo -e "Non-optional parameters for fasta files or reads were not entered\n"
#	echo "Attempting using the expected automatic parameters assumed by MuDoGeR"
#	
	
#fi

#Checking if all non-optional parameters are entered correctly
if [ "$reduced" = "false" ] && [ "$complete" = "false" ]; then 
	echo "Non-optional parameter for BRAT type was not entered. Please select --reduced or --complete."
	help_message; exit 1
fi

#Checking if all non-optional parameters are entered correctly
if [ "$reduced" = "true" ] && [ "$complete" = "true" ]; then 
	echo "Non-optional parameter for BRAT type was not entered. Please only select one option, --reduced or --complete."
	help_message; exit 1
fi

#Checking if all non-optional parameters are entered correctly
if [ "$absolute" = "false" ] && [ "$coverage" = "false" ] && [ "$relative" = "false" ]; then 
	echo "Non-optional parameters for output type were not entered. Please only select one option, --absolute-values ,--coverage or--relative-abundance."
	help_message; exit 1
fi

#Checking if output directory already exists
if [ ! -d $WORKDIR/mapping_results ]; then mkdir -p $WORKDIR/mapping_results; fi


echo -e "\nTHE PIPELINE STARTED\n"


#Module 5 - BRATs for MuDoGer

# FASTAS = bin folder
# READS = lib folder
# WORKDIR = output folder
# THREAD =number of threads


conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path

#     OTUpick                         submodule 5-1
#     Mapping process	              submodule 5-2

if [ "$brat_type" = "--complete" ]; then

	echo -e "\n EUK BIN CALCULATION STARTED"

	mkdir -p "$libname_folder"/eukaryotes/
	bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-4-1_eukrep-eukbin-filter.sh "$assembly"       \
                                      "$forward_library"                                 \
                                      "$reverse_library"                                 \
                                      "$libname_folder"/eukaryotes			 \
                                      "$cores"                                           \  
				      "$memory"

	echo -e "\n EUK BIN CALCULATION DONE"
	
elif [ "$brat_type" = "--reduced" ]; then

	#     OTU picking      submodule 5-1
	echo -e "\n OTU picking STARTED"

	#if [ -z "$(ls -A "$libname_folder"/eukaryotes/filtered_euk_bins/)" ]; then
   	#	echo -e "\nNo relevant eukaryotic found"; touch "$libname_folder"/eukaryotes/no_euk_bins_for_genemark
	#else
   		bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-5-1_prokOTUpicking.sh $WORKDIR $metadata_table $THREAD
	#fi

	echo -e "\n OTU picking DONE"
	
fi




exit 1










#1.0 index bins

#echo "                                      STARTING: 1.0 index bins"

#cd $FASTA
#for bin in *.fa ; do echo "${bin/.fa/}"; done > aux
#parallel -j $THREAD /data/msb/tools/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2-build {}.fa {} < aux
#rm -f aux

# old command
#for bin in *.fa ; do /data/msb/tools/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2-build "$bin" "${bin/.fa/}"; done

# rm later all bins with *.bt2


#1.1 calculate genome size for coverage
if [ "$coverage" = "true" ]; then

echo "                                      STARTING: 1.1 calculate genome size for coverage"

mkdir temp.genome-size

for bin in *.fa; do grep -v ">" "$bin" > aux; size="$(wc -c aux | cut -f1 -d' ' )"; echo $size "${bin/.fa/}" > temp.genome-size/"${bin/.fa/}.genome-size"; rm -f aux; done
 
cd "$FASTA"/temp.genome-size
cat *.genome-size > "$WORKDIR"/genomes-size.txt
#rm -fr "$FASTA"/temp.genome-size

fi

#(deactivated ) 2.0 merge PE reads from each library to fasta file
# deactivated 
# cd $READS
# mkdir merged

# module load GCCcore/6.4.0 zlib/1.2.11 bzip2/1.0.6 libtool/2.4.6

# for d in *_1.fastq ; do forward=$READS/"$d"; reverse="${forward/_1./_2.}";  /gpfs1/data/msb/tools/pandaseq/pandaseq/pandaseq -f "$forward" -r "$reverse" -w "$READS"/merged/"${d/_1.fastq/.fasta}"

# convert SE reads to fasta file. All file need to be *_SE.fastq

# for SE_lib in *_SE.fastq; do sed -n '1~4s/^@/>/p;2~4p' $SE_lib > merged/"${SE_lib/_SE.fastq/.fasta}"

# 2.1) Count number of reads in the libs (removed

##### The number of hits is already counted in the “frag size”. So remove this step was removed.

#mkdir $READS/total_reads
#mkdir -p $READS/merged/total_reads
#module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.0 Python/3.7.4

# Count for PE reads. 
# for i in $READS/*_1.fastq; do lib="$(echo $i | rev | cut -f1 -d "/" | rev | sed "s/_1.fastq//g")" ; wc -l $i > "$READS"/total_reads/"$lib".test1 ; awk '{print $2, $1/4}' "$READS"/"$lib".test1  | sed "s/_1.fastq//g" > "$READS"/total_reads/"$lib".total-reads; rm -f "$READS"/"$lib".test1

#Count for SE reads
#for i in $READS/*_SE.fastq; do lib="$(echo $i | rev | cut -f1 -d "/" | rev | sed "s/_SE.fastq//g")" ; wc -l $i > "$READS"/total_reads/"$lib".test1 ; awk '{print $2, $1/4}' "$READS"/"$lib".test1  | sed "s/_SE.fastq//g" > "$READS"/total_reads/"$lib".total-reads; rm -f "$READS"/"$lib".test1

# count for merged reads in fasta file. Number of lines divided by 2 as there is not a quality score lines
#for i in $READS/merged/*.fasta; do lib="$(echo $i | rev | cut -f1 -d "/" | rev | sed | sed "s/.fasta//g")" ; wc -l $i > "$READS"/merged/total_reads/"$lib".test1 ; awk '{print $2, $1/2}' "$READS"/merged/"$lib".test1 | sed "s/.fasta//g" > "$READS"/merged/total_reads/"$lib".total-reads ; rm -f  "$READS"/merged/"$lib".test1 
# Combine results. 
#mkdir "$WORKDIR"/total_reads
#cat "$READS"/merged/total_reads/*.total-reads > "$WORKDIR"/total_reads/total_reads_per_lib.txt
#rm -fr "$READS"/merged/total_reads
#cat "$READS"/total_reads/* >> "$WORKDIR"/total_reads/total_reads_per_lib.txt 
#rm -fr "$READS"/total_reads
#2.1) Frag size and total number of reads of the libs for coverage and relative abundance

if [ "$coverage" = "true" ] || [ "$relative" = "true"  ]; then

# get fragment size. First remove the headers, then it is the size of the lib (number of characters (wc -c)) divided by the number of sequences (number of lines(wc -l)).
# for PE reads

echo "                                      STARTING: Frag size and total number of reads of the libs for coverage and relative abundance"

cd "$READS"
mkdir "$WORKDIR"/frag-size "$WORKDIR"/total_reads 

# How to convert from fastq to fasta
# for lib in *_1.fastq; do sed -n '1~4s/^@/>/p;2~4p' $lib > temp/"${lib/_1.fastq/_1.fasta}"; done
# First transform the forward reads into a sequence only of the reads, no headers or quality info 
#convert to parallelized jobs
mkdir  temp
# the command below to get only the reads (line 2 of each 4 lines)
for lib in *_1.fastq; do echo " sed -n '2~4p' $lib > "$WORKDIR"/temp/${lib/_1.fastq/_1.only-reads} "; done > to_only_reads.txt

parallel -j $THREAD {} < to_only_reads.txt

#rm -f to_only_reads.txt

# parallel -j 2 {} < to_only_reads.txt

cd "$WORKDIR"/temp
# count number of reads and of frag-size
## for lib in *.only-reads; do seqs_num="$(wc -l $lib | cut -f1 -d' ' )" ; size="$(wc -c $lib | cut -f1 -d' ' )"; frag_avg_size="$((size/seqs_num))"; echo $frag_avg_size "${lib/_1.only-reads/}" > "$WORKDIR"/frag-size/"${lib/_1.only-reads/}.frag-size"; echo "$seqs_num" "${lib/_1.only-reads/}" > "$WORKDIR"/total_reads/"${lib/_1.only-reads/}.total_reads" ; done 
# it is not working properly because it is calculating before printing the commands, so I am running one by one for now
# for lib in *.only-reads; do echo "seqs_num=\"$(wc -l $lib | cut -f1 -d' ' )\" ; size=\"$(wc -c $lib | cut -f1 -d' ' )\"; frag_avg_size=\"$((size/seqs_num))\"; echo $frag_avg_size \"${lib/_1.only-reads/}\" > \"$WORKDIR\"/frag-size/\"${lib/_1.only-reads/}.frag-size\"; echo \"$seqs_num\" \"${lib/_1.only-reads/}\" > \"$WORKDIR\"/total_reads/\"${lib/_1.only-reads/}.total_reads\" " ; done > to_frag_size_and_total_reads.txt

#parallel -j $THREAD {} < to_frag_size_and_total_reads.txt

# parallel -j 1 {} < to_frag_size_and_total_reads.txt

for lib in *.only-reads; do seqs_num="$(wc -l $lib | cut -f1 -d' ' )" ; size="$(wc -c $lib | cut -f1 -d' ' )"; frag_avg_size="$((size/seqs_num))"; echo $frag_avg_size "${lib/_1.only-reads/}" > "$WORKDIR"/frag-size/"${lib/_1.only-reads}".frag-size; echo "$seqs_num" "${lib/_1.only-reads}" > "$WORKDIR"/total_reads/"${lib/_1.only-reads}".total_reads; done 

# old commands
# for lib in *_1.fastq; do grep -A1 ">" "$lib" |  grep -v ">" > aux; seqs_num="$(wc -l aux | cut -f1 -d' ' )" ; size="$(wc -c aux | cut -f1 -d' ' )"; frag_avg_size="$((size/seqs_num))"; echo $frag_avg_size "${lib/_1.fastq/}" > "$WORKDIR"/frag-size/"${lib/_1.fastq/}.frag-size"; echo "$seqs_num" "${lib/_1.fastq/}" > "$WORKDIR"/total_reads/"${lib/_1.fastq/}.total_reads" ; rm -f aux; done 

# for lib in *_1.fastq; do grep -A1 ">" "$lib" |  grep -v ">" > aux; seqs_num="$(wc -l aux | cut -f1 -d' ' )" ; size="$(wc -c aux | cut -f1 -d' ' )"; frag_avg_size="$((size/seqs_num))"; echo $frag_avg_size "${lib/_1.fastq/}" > "$WORKDIR"/frag-size/"${lib/_1.fastq/}.frag-size"; echo "$seqs_num" "${lib/_1.fastq/}" > "$WORKDIR"/total_reads/"${lib/_1.fastq/}.total_reads" ; rm -f aux; done 
# testing parallelization
#mkdir temp
# First, change files to a new format only with reads and no header.
# for lib in *_1.fastq; do echo "grep -A1 \">\" "$lib" |  grep -v \">\"  > temp/"$lib".aux" ; done > to_fasta.txt

# parallel -j $THREAD {} < to_fasta.txt
# for lib in *_1.fastq; do echo grep -A1 \">\" "$lib" |  grep -v \">\"  > temp/"$lib".aux ; done > 
# for lib in *_1.fastq; do grep -A1 ">" "$lib" |  grep -v ">" ; done
# seqs_num="$(wc -l aux | cut -f1 -d' ' )" ; size="$(wc -c aux | cut -f1 -d' ' )"; frag_avg_size="$((size/seqs_num))"; echo $frag_avg_size "${lib/_1.fastq/}" > "$WORKDIR"/frag-size/"${lib/_1.fastq/}.frag-size"; echo "$seqs_num" "${lib/_1.fastq/}" > "$WORKDIR"/total_reads/"${lib/_1.fastq/}.total_reads" ; rm -f aux; done 
# for lib in *_1.fasta; do grep -v ">" "$lib" > aux; seqs_num="$(wc -l aux | cut -f1 -d' ' )" ; size="$(wc -c aux | cut -f1 -d' ' )"; frag_avg_size="$((size/seqs_num))"; echo $frag_avg_size "${lib/_1.fastq/}" > "$WORKDIR"/frag-size/"${lib/_1.fastq/}.frag-size"; echo "$((seqs_num*2))" "${lib/_1.fastq/}" > "$WORKDIR"/total_reads/"${lib/_1.fastq/}.total_reads" ; rm -f aux; done 
# Combine results

cd "$WORKDIR"/total_reads/
cat *.total_reads > "$WORKDIR"/total_reads_per_lib.txt 

cd "$WORKDIR"/frag-size/
cat *.frag-size > "$WORKDIR"/frag-size_per_lib.txt 

# cd "$WORKDIR"

#rm -fr "$WORKDIR"/total_reads "$WORKDIR"/frag-size

# for lib in *_1.fasta; do grep -v ">" "$lib" > aux; seqs_num="$(wc -l aux | cut -f1 -d' ' )" ; size="$(wc -c aux | cut -f1 -d' ' )"; frag_avg_size="$((size/seqs_num))"; echo $frag_avg_size "${lib/_1.fastq/}" > frag-size/"${lib/_1.fastq/}.frag-size" ; rm -f aux; done 

#cd "$READS"/frag-size
#cat *.frag-size > $WORKDIR/lib-frag-sizes

# for SE reads
#cd "$READS"
#for lib in *_SE.fasta; do grep -v ">" "$lib" > aux; seqs_num="$(wc -l aux | cut -f1 -d' ' )" ; size="$(wc -c aux | cut -f1 -d' ' )"; frag_avg_size="$((size/seqs_num))"; echo $frag_avg_size "${lib/_SE.fastq/}" > frag-size/"${lib/_SE.fastq/}.frag-size"; rm -f aux; done 
#cd "$READS"/frag-size
#cat *.frag-size > $WORKDIR/lib-frag-sizes

# for reads in fasta
#cd "$READS"/merged/
#for lib in *.fasta; do grep -v ">" "$lib" > aux; seqs_num="$(wc -l aux | cut -f1 -d' ' )" ; size="$(wc -c aux | cut -f1 -d' ' )"; frag_avg_size="$((size/seqs_num))"; echo $frag_avg_size "${lib/.fasta/}" > "$READS"/frag-size/"${lib/.fasta/}.frag-size"; rm -f aux; done 
#cd $WORKDIR/libs/merged/frag-size
#cat *.frag-size > $WORKDIR/lib-frag-sizes
fi


# 3) map the libraries to the indexed bins
mkdir "$WORKDIR"/temp
cd $FASTA

if [ "$complete" = "true" ]; then

echo "                                      STARTING: map the libraries to the indexed bins '--complete'"
module load foss/2020b metaWRAP/1.3-Python-2.7.18

# Complete BRATs input job list
for d in *.fa; do bin="$FASTA"/"${d/.fa/}"; for l in "$READS"/*_1.fastq; do lib="${l/_1.fastq/}";  output="$WORKDIR"/temp/"${d/.fa/}"-LIB-"$(echo $lib | rev | cut -f1 -d'/' | rev)".txt; echo "$bin" "$lib"_1.fastq "$lib"_2.fastq "$output" ; done; done > "$WORKDIR"/mapping_job_list_PE-complete.txt

# Create mapping complete BRATs/ mapping job list


cat "$WORKDIR"/mapping_job_list_PE-complete.txt | while read l ; 
do bin="$(echo $l | cut -f1 -d " " )" ; lib1="$(echo $l | cut -f2 -d " " )" ; lib2="$(echo $l | cut -f3 -d " " )" ; out="$(echo $l | cut -f4 -d " " )" ; echo "/gpfs1/data/msb/tools/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2 --threads 1 -x" "$bin" -1 "$lib1" -2 "$lib2" "|" "samtools view  --threads 1 -F 4 |  cut -f1 |  sort -u | wc -l" ">" "$out" ; done > aux1

# run mapping jobs 

parallel -j $THREAD {} < aux1

#rm -f aux1 "$WORKDIR"/mapping_job_list_PE-complete.txt

fi
if [ "$reduced" = "true" ]; then

echo "                                      STARTING: map the libraries to the indexed bins '--reduced'"
module load foss/2020b metaWRAP/1.3-Python-2.7.18

# Reduced BRATs input job list
for d in *.fa; do bin="$FASTA"/"${d/.fa/}"; lib="$(echo $d | sed "s/-bin./\t/g" | cut -f1 )" ; output="$WORKDIR"/temp/"${d/.fa/}"-LIB-"$(echo $lib | rev | cut -f1 -d'/' | rev)".txt; echo "$bin" "$lib"_1.fastq "$lib"_2.fastq "$output" ; done > "$WORKDIR"/mapping_job_list_PE-reduced.txt
# for d in *.fa; do bin="$FASTA"/"${d/.fa/}"; for l in "$READS"/*_1.fastq; do lib="${l/_1.fastq/}";  output="$WORKDIR"/temp/"${d/.fa/}"-LIB-"$(echo $lib | rev | cut -f1 -d'/' | rev)".txt; echo "$bin" "$lib"_1.fastq "$lib"_2.fastq "$output" ; done; done > "$WORKDIR"/mapping_job_list_PE.txt

# Create SE reads job array 
# for d in *.fa; do bin="$FASTA"/"${d/.fa/}"; for l in "$READS"/*_SE.fastq; do lib="${l/_SE.fastq/}";  output="$WORKDIR"/temp/"${d/.fa/}"-LIB-"$(echo $lib | rev | cut -f1 -d'/' | rev)".txt; echo "$bin" "$lib"_SE.fastq "$output" ; done; done > "$WORKDIR"/mapping_job_list_SE.txt

# Create merged reads job array 
# for d in *.fa; do bin="$FASTA"/"${d/.fa/}"; for l in "$READS"/merged/*.fasta; do lib="${l/.fastq/}";  output="$WORKDIR"/temp/"${d/.fa/}"-LIB-"$(echo $lib | rev | cut -f1 -d'/' | rev)".txt; echo "$bin" "$lib".fastq "$output" ; done; done > "$WORKDIR"/mapping_job_list_merged.txt

# Create mapping reduced BRATs/ mapping job list
# I need to change to the new reduced using the groups from gOTUpick.

cat "$WORKDIR"/mapping_job_list_PE-reduced.txt | while read l ; 
do bin="$(echo $l | cut -f1 -d " " )" ; lib1="$(echo $l | cut -f2 -d " " )" ; lib2="$(echo $l | cut -f3 -d " " )" ; out="$(echo $l | cut -f4 -d " " )" ; echo "/gpfs1/data/msb/tools/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2 --threads 1 -x" "$bin" -1 "$lib1" -2 "$lib2" "|" "samtools view  --threads 1 -F 4 |  cut -f1 |  sort -u | wc -l" ">" "$out" ; done > aux2

# run mapping jobs. aux1 is complete BRAT and aux2 is reduced BRAT 

parallel -j $THREAD {} < aux2

#rm -f aux2


#cat "$WORKDIR"/mapping_job_list_PE.txt | while read l ; 
#do bin="$(echo $l | cut -f1 -d " " )" ; 
#lib1="$(echo $l | cut -f2 -d " " )" ; 
#lib2="$(echo $l | cut -f3 -d " " )" ; 
#out="$(echo $l | cut -f4 -d " " )" ;
#/gpfs1/data/msb/tools/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2 \
#        --threads ${SLURM_CPUS_PER_TASK:-1} \
#        -x "$bin" \
#        -1 "$lib1" \
#        -2 "$lib2" | samtools view  --threads ${SLURM_CPUS_PER_TASK:-1} -F 4 |
#        cut -f1 |
#        sort -u |
#        wc -l > "$out" ;
#	done


# Count hits with SE reads 
#cat "$WORKDIR"/mapping_job_list_SE.txt | while read l ; 
#do bin="$(echo $l | cut -f1 -d " " )" ; 
#lib="$(echo $l | cut -f2 -d " " )" ;
#out="$(echo $l | cut -f3 -d " " )";
#/gpfs1/data/msb/tools/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2 \
#        --threads ${SLURM_CPUS_PER_TASK:-1} \
#        -x "$bin" \
#        -q "$lib" | samtools view --threads ${SLURM_CPUS_PER_TASK:-1}  -F 4 |
#        cut -f1 |
#        sort -u |
#        wc -l > "$out" ;
#	done

# Count hits with fasta file
#cat "$WORKDIR"/mapping_job_list_merged.txt | while read l ; 
#do bin="$(echo $l | cut -f1 -d " " )" ; 
#lib="$(echo $l | cut -f2 -d " " )" ;
#out="$(echo $l | cut -f3 -d " " )"; /gpfs1/data/msb/tools/bowtie2/bowtie2-2.4.1-linux-x86_64/bowtie2 \
#        --threads ${SLURM_CPUS_PER_TASK:-1} \
#        -x "$bin" \
#        -f "$lib" | samtools view --threads ${SLURM_CPUS_PER_TASK:-1}  -F 4 |
#        cut -f1 |
#        sort -u |
#        wc -l > "$out" ;
#	done

fi

# 4) BRATs absolute number of mapping
# merge results. Edit the sed according to the parameters that separate the bin and the lib

echo "                                      STARTING: 4) BRATs number of mapping"

cd $WORKDIR/temp
for d in *.txt; do echo -e "$d\t\c" | sed "s/-LIB-/\t/g" | sed "s/.txt//g"  ; cat $d; done > "$WORKDIR"/BRAT_absolute_values_list.tsv

# create cross-table (I had problem with that in one project but it is what we have to run in bash)
cat "$WORKDIR"/BRAT_absolute_values_list.tsv | /data/msb/tools/datamash/datamash-1.3/datamash -sW crosstab 1,2 unique 3 > "$WORKDIR"/BRAT_absolute_values_cross-table.tsv


# 4.1) Coverage


if [ "$coverage" = "true" ]; then

echo "                                      STARTING: 4.1) Coverage"

#"$WORKDIR"/total_reads_per_lib.txt 
#"$WORKDIR"/frag-size_per_lib.txt 

#"$WORKDIR"/genomes-size.txt
#"$WORKDIR"/BRAT_absolute_values_list.tsv
# I need to check the separators for the cuts.

cd "$WORKDIR"
while read l; do num_hits="$(echo $l | cut -f3 -d " ")" ; lib="$(echo $l | cut -f2 -d " ")"  ; bin="$(echo $l | cut -f1 -d " ")" ; frag_size="$(grep -w $lib frag-size_per_lib.txt | cut -f2 )" ; gen_size="$(grep -w $bin genomes-size.txt | cut -f1 -d ' ')" ; hits_times_frag="$( bc -l <<< "$num_hits*$frag_size" )"; coverage="$( bc -l <<< "$hits_times_frag/$gen_size" )"; echo -e "$bin""\t""$lib""\t""$coverage"; done < BRAT_absolute_values_list.tsv > BRAT_coverage_list.tsv

# create cross-table (I had problem with that in one project but it is what we have)
cat "$WORKDIR"/BRAT_coverage_list.tsv | /data/msb/tools/datamash/datamash-1.3/datamash -sW crosstab 1,2 unique 3 > "$WORKDIR"/BRAT_coverage_cross-table.tsv

if [ ! -s "$WORKDIR"/BRAT_coverage_cross-table.tsv ]; then echo "Something went wrong while creating the final output file. Exiting..."; exit 1; fi 
# while read l; do num_hits="$(echo $l | cut -f3 -d " ")" ; lib="$(echo $l | cut -f2 -d " ")"  ; bin="$(echo $l | cut -f1 -d " ")" ; frag_size="$(grep -w $lib frag-size_per_lib.txt | cut -f2 )" ; gen_size="$(grep -w $bin genomes-size.txt | cut -f1 -d ' ')" ; hits_times_frag="$(( $num_hits * $frag_size))"; coverage="$((hits_times_frag/gen_size))"; echo "$bin" "$lib" "$coverage"; done < BRAT_absolute_values_list.tsv > BRAT_coverage_list.tsv

# while read l; do num_hits="$(echo $l | cut -f3 -d " ")" ; lib="$(echo $l | cut -f2 -d " ")"  ; bin="$(echo $l | cut -f1 -d " ")" ; frag_size="$(grep -w $lib frag-size_per_lib.txt | cut -f2  )" ; gen_size="$(grep -w $bin genomes-size.txt | cut -f1 -d ' ')" ; echo $frag_size; done < BRAT_absolute_values_list.tsv
fi
# 4.2) Relative Abundance
if [ "$relative" = "true" ]; then

echo "                                      STARTING: 4.2) Relative Abundance"

while read l; do num_hits="$(echo $l | cut -f3 -d " ")" ; lib="$(echo $l | cut -f2 -d " ")"  ; bin="$(echo $l | cut -f1 -d " ")" ; n_reads="$(grep -w $lib total_reads_per_lib.txt | cut -f2 -d " " )" ; relative_abundance="$( bc -l <<< "$num_hits/$n_reads" )"; echo -e "$bin""\t""$lib""\t""$relative_abundance"; done < BRAT_absolute_values_list.tsv > BRAT_relative_abundance_list.tsv

cat "$WORKDIR"/BRAT_relative_abundance_list.tsv | /data/msb/tools/datamash/datamash-1.3/datamash -sW crosstab 1,2 unique 3 > "$WORKDIR"/BRAT_relative_abundance_cross-table.tsv

if [ ! -s "$WORKDIR"/BRAT_relative_abundance_cross-table.tsv ]; then echo "Something went wrong while creating the final output file. Exiting..."; exit 1; fi 

fi

#if [ "$absolute" = "false" ]; then

#rm -f "$WORKDIR"/BRAT_absolute_values_list.tsv "$WORKDIR"/BRAT_absolute_values_cross-table.tsv

#fi

echo -e "\nThat's the end of BRAT. Enjoy your results. Goodbye!"
