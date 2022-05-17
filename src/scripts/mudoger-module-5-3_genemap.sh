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
mkdir -p $genelength_results_path
mkdir -p $genemap_count_results_path


# loop around samples
cd $genemap_results_path

aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
for i in $aux; 
do

#Run bowtie2-build
bowtie2-build $WORKDIR/$i/$assembly_input_file $i.reference;

#Run mapping with bowtie2
bowtie2 -p $cores -x $i.reference -1 $WORKDIR/$i/$qc_input_path/final_pure_reads_1.fastq -2 $WORKDIR/$i/$qc_input_path/final_pure_reads_2.fastq -S $i.map.sam;

#Convert SAM to BAM and sort
samtools sort -o $i.map.sorted.bam -@ $cores -O bam $i.map.sam;

done

cd -
conda deactivate

#Run functional annotation with prokka
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/prokka_env

aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
for i in $aux; 
do

#Run prokka on assembly
prokka $WORKDIR/$i/$assembly_input_file --outdir $functional_assembly_path/$i --prefix $i --metagenome --cpus $cores;

#Covert gff to gtf
grep -v "#" $functional_assembly_path/$i/$i.gff | grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,4,5,7,9 |  awk -v OFS='\t' '{print $1,"PROKKA","CDS",$2,$3,".",$4,".","gene_id " $5}' > $functional_assembly_path/$i/$i.gtf

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

### Missing covarege and TPM calculation ####

