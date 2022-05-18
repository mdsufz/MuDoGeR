#!/bin/bash

#Steps
#Assembly bowtie index and qc reads mapping
#Assembly gene annotation by prokka
#Convert .gff to .gtf file
#Count absolute number of reads mapped

#The user should have initial data to use this script: Assembly fasta file, Quality-controlled reads.

conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

#Input parameters
WORKDIR="$1"
metadata_table="$2"
cores="$3"
absolute=$4;
coverage=$5;
relative=$6;

#Define dependent parameters
assembly_input_file="assembly/final_assembly.fasta"
qc_input_path="qc/"

#Define output files path
genemap_results_path="$WORKDIR/mapping_results/assembly_gene_map/raw_mapping"
functional_assembly_path="$WORKDIR/mapping_results/assembly_gene_map/functional_annotation/"
genelength_results_path="$WORKDIR/mapping_results/assembly_gene_map/genelength/"
genemap_count_results_path="$WORKDIR/mapping_results/assembly_gene_map/map_absolute_count/"
genemap_cov_results_path="$WORKDIR/mapping_results/assembly_gene_map/map_coverage_norm/"
genemap_tpm_results_path="$WORKDIR/mapping_results/assembly_gene_map/map_tpm_norm/"

#Gene mapping started
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/brat_env

mkdir -p $genemap_results_path
mkdir -p $functional_assembly_path
mkdir -p $genemap_count_results_path

# loop around samples
cd $genemap_results_path

aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
for i in $aux; 
do

#Run bowtie2-build

if [ -f  $genemap_results_path/$i.reference.1.bt2 ];
then echo "-> Bowtie-build already done. Please check here: "$genemap_results_path"/$i.reference"
else
echo "-> Running bowtie-build"

bowtie2-build $WORKDIR/$i/$assembly_input_file $i.reference;

fi

#Run mapping with bowtie2
if [ -f  $genemap_results_path/$i.map.sam ];
then echo "-> Map already done. Please check here: $genemap_results_path/$i.map.sam"
else
echo "-> Mapping reads started"

bowtie2 -p $cores -x $i.reference -1 $WORKDIR/$i/$qc_input_path/final_pure_reads_1.fastq -2 $WORKDIR/$i/$qc_input_path/final_pure_reads_2.fastq -S $i.map.sam;

fi

#Convert SAM to BAM and sort
if [ -f  $genemap_results_path/$i.map.sorted.bam ];
then echo "-> Sorting .sam file already done. Please check here: $genemap_results_path/$i.map.sorted.bam"
else
echo "-> Sorting .sam file"

samtools sort -o $i.map.sorted.bam -@ $cores -O bam $i.map.sam;

fi

done

cd -
conda deactivate

#Run functional annotation with prokka
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/prokka_env

aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
for i in $aux; 
do


#Run prokka on assembly
if [ -f  $functional_assembly_path/$i/$i.gff ];
then echo "-> Gene annotation already done. Please check here: $functional_assembly_path/$i/$i.gff"
else
echo "-> Gene annotation started"

prokka $WORKDIR/$i/$assembly_input_file --outdir $functional_assembly_path/$i --prefix $i --metagenome --cpus $cores;

fi

#Covert .gff to .gtf
if [ -f  $functional_assembly_path/$i/$i.gtf ];
then echo "-> GTF file already done. Please check here: $functional_assembly_path/$i/$i.gtf"
else
echo "-> Converting to GTF file"

grep -v "#" $functional_assembly_path/$i/$i.gff | grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,4,5,7,9 |  awk -v OFS='\t' '{print $1,"PROKKA","CDS",$2,$3,".",$4,".","gene_id " $5}' > $functional_assembly_path/$i/$i.gtf

fi

done

conda deactivate

#Run count the reads mapped
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/htseq_env

cd $genemap_results_path

aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
for i in $aux; 
do

if [[ -f "$genemap_count_results_path/$i.count" ]]; then
echo "$genemap_count_results_path/$i.count already exists"
continue
else
echo "Counting $i started"
htseq-count -r pos -t CDS -f bam $i.map.sorted.bam $functional_assembly_path/$i/$i.gtf > $genemap_count_results_path/$i.count
fi

done
cd -

conda deactivate

#Calculating Coverage

if [ "$coverage" = "true" ] || [ "$relative" = "true" ]; then

  conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/cov_env

  mkdir -p $genelength_results_path
  mkdir -p $genemap_cov_results_path
  if [ "$relative" = "true" ]; then
    mkdir -p $genemap_tpm_results_path
  fi

  cd $genemap_cov_results_path

  rm -f $WORKDIR/mapping_results/assembly_gene_map/avg_reads_len.tsv
  aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
  for i in $aux; 
  do

  #Calculate average read length from samples
  echo "Calculating average read lenght from $i"
  avg_len=`awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' $WORKDIR/$i/$qc_input_path/final_pure_reads_1.fastq`
  echo -e "$i\t$avg_len" >> $WORKDIR/mapping_results/assembly_gene_map/avg_reads_len.tsv

  #Calculate gene length
  if [[ -f "$genelength_results_path/$i.genelength" ]]; then
  echo "$genelength_results_path/$i.genelength already exists"
  continue
  else

  cut -f4,5,9 $functional_assembly_path/$i/$i.gtf | sed 's/gene_id //g' | gawk '{print $3,$2-$1+1}' | tr ' ' '\t' > $genelength_results_path/$i.genelength

  fi

  #Calculate gene coverage and TPM in the sample
  if [[ -f "$genemap_cov_results_path/$i.cov" ]]; then
  echo "$genemap_cov_results_path/$i.cov already exists"
  continue
  else

  avg_len=$(cat $WORKDIR/mapping_results/assembly_gene_map/avg_reads_len.tsv | grep "$i"); python "$MUDOGER_DEPENDENCIES_PATH"/tpm_cov_table_mudoger.py -n $i -c $genemap_count_results_path/$i.count -i <(echo -e "$avg_len") -l $genelength_results_path/$i.genelength

  mv $i.cov $genemap_cov_results_path/$i.cov
  if [ "$relative" = "true" ]; then
    mv $i.tpm $genemap_tpm_results_path/$i.tpm
  elif [ "$relative" = "false" ]; then
    rm -f $i.tpm
  fi
  
  fi

  done

  cd -
  conda deactivate
fi
