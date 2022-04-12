#!/bin/bash

# this bash script generates some simple genome statistics and summary from Uvigs analysis
# should be used like: N50.sh Multi_fasta_file

##### Base of the script: https://github.com/hcdenbakker/N50.sh/blob/master/N50.sh
# loading conda environment
echo '------- START MODULE 3-5 Uvigs METRICS'
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database


viruses_folder="$1"viruses
output_folder="$1"


#path_to_lib='/Users/toscan/Documents/trabalho/ufz/mudoger/A/'

viruses_folder=$path_to_lib'viruses'

#lib="$(echo $path_to_lib | rev | cut -f2 -d'/' | rev )"
host_results=$viruses_folder'/host_prediction/output_results'
derep=$viruses_folder'/investigation/dereplication'
uvigs=$viruses_folder'/host_prediction/uvigs'
quality_summary=$viruses_folder'/vcheck_quality/quality_summary.tsv'

if [ ! -f $derep'/uvigs_mapping.txt' ];
then for uvig in $uvigs/*; do echo -e "$(echo $uvig | rev | cut -f1 -d'/' | rev )\t\c"; echo "$(cat $uvig | grep '>' | sed "s/>//g" )" ; done > $derep'/uvigs_mapping.txt'
else :; fi

if [ !  -f $viruses_folder'/viruses_summary.tsv' ]; 
then
while read l; do
uvig="$(echo "$l" | cut -f1 | cut -f1 -d'.')"
contig="$(echo "$l" | cut -f2)"
echo -e "$uvig\t\c"
echo -e "$(grep $contig $quality_summary)\t\c"
grep -w $uvig $host_results'/prediction.list' | cut -f2,3
done < $derep'/uvigs_mapping.txt' > $viruses_folder'/.viruses_summary_raw.tsv'
echo -e 'uvig\toriginal_contig\tuvig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings\tputative_host\tlikelihood' >  $viruses_folder'/.header'
cat $viruses_folder'/.header' $viruses_folder'/.viruses_summary_raw.tsv' > $viruses_folder'/viruses_summary.tsv'
else :; 
fi

# FIlter Good quality bins (Complet - 5*Contamination >= 50) based on CheckM
#cat $output_path/genome_metrics.tsv | awk 'BEGIN {FS="\t"};NR==1; ($2 - (5*$3)) >= 50' > "$1"/MAGS_results.tsv


