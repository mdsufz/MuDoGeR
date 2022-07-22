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


viruses_folder="$1"/viruses
output_folder="$1"


host_results=$viruses_folder'/host_prediction/output_results'
derep=$viruses_folder'/investigation/dereplication'
uvigs=$viruses_folder'/final_outputs/putative_viral_contigs'
quality_summary=$viruses_folder'/vcheck_quality/quality_summary.tsv'
vcontact_output=$viruses_folder'/taxonomy/vcontact-output'

mkdir -p "$viruses_folder"/final_outputs/

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

rm -f $viruses_folder'/.header'
rm -f $viruses_folder'/.viruses_summary_raw.tsv'

echo "Formating Vcontact2 results..."

if [ !  -f $viruses_folder'/taxa_estimation_vir_vcontact2_summary.csv' ];
then
    	  rm -f $viruses_folder'/.taxa_viruses_summary_raw.csv'
        while read l;
          do
          uvig="$(echo "$l" | cut -f1 | cut -f1 -d'.')"
          contig="$(echo "$l" | cut -f2)"
          contig_line=$(echo | cat $vcontact_output/genome_by_genome_overview.csv | grep "$contig")
          clt="$(echo "$contig_line" | cut -f7 -d',')"
          if [ "$clt" == "Clustered" ];
                then
                vc_clt=$(echo "$contig_line" | cut -f9 -d',')
                cat $vcontact_output/genome_by_genome_overview.csv | grep "$vc_clt" | cut -f3,4 -d',' > aux
                while read j;
                  do
                  ord="$(echo $j | cut -f1 -d',')";
                  if [ "$ord" != "Unassigned" ];
                    then
                    ord_family=$j
                    break;
                  else
                    ord_family="Virus,Unassigned,Unassigned";
                  fi;
                done < aux
                rm -f aux
          else
               	ord_family="Virus,Unassigned,Unassigned"
          fi
          echo -e "$uvig,$contig,$ord_family" >> $viruses_folder'/.taxa_viruses_summary_raw.csv'
        done < $derep'/uvigs_mapping.txt'
echo -e 'uvig,original_contig,Domain,Order,Family' >  $viruses_folder'/.taxa_header'
cat $viruses_folder'/.taxa_header' $viruses_folder'/.taxa_viruses_summary_raw.csv' > $viruses_folder'/taxa_estimation_vir_vcontact2_summary.csv'
else :; 
fi
rm -f $viruses_folder'/.taxa_header'
rm -f $viruses_folder'/.taxa_viruses_summary_raw.csv'

#Filter Good quality Uvigs based on CheckV
awk -F "\t" 'NR==1;{ if($10 == "High-quality") { print } }' $viruses_folder'/viruses_summary.tsv' > $viruses_folder/uvigs_high_quality.tsv

# Organize results in a single results folder

mkdir -p "$viruses_folder"/final_outputs/
mkdir -p "$viruses_folder"/final_outputs/only_uvigs_seq
mkdir -p "$viruses_folder"/final_outputs/putative_vir_seq_metrics_summary # Hosts, Quality, Taxa, Genes

#Move uvigs summary file
mv -f $viruses_folder/uvigs_high_quality.tsv "$viruses_folder"/final_outputs/uvigs_high_quality_summary.tsv
#Move summary file
mv -f $viruses_folder'/viruses_summary.tsv' "$viruses_folder"/final_outputs/putative_vir_contigs_summary.tsv
#Move mapping file
#mv $derep'/uvigs_mapping.txt' "$viruses_folder"/final_outputs/viral_contigs_seq_names.txt
sed "s/\t/,/g" $derep'/uvigs_mapping.txt' > "$viruses_folder"/final_outputs/viral_contigs_seq_names.csv

#Copy only UViGs seq
cat "$viruses_folder"/final_outputs/uvigs_high_quality_summary.tsv | cut -f1 | tail -n+2 > "$viruses_folder"/final_outputs/tmp
while read uvig; do yes | cp "$viruses_folder"/final_outputs/putative_viral_contigs/$uvig.fa "$viruses_folder"/final_outputs/only_uvigs_seq;
done < "$viruses_folder"/final_outputs/tmp

rm -f "$viruses_folder"/final_outputs/tmp

#Copy Host list file summary
yes | cp $host_results'/prediction.list' "$viruses_folder"/final_outputs/putative_vir_seq_metrics_summary/host_vir_pair_wish_summary.tsv

#Copy quality file summary
cat "$viruses_folder"/final_outputs/putative_vir_contigs_summary.tsv | cut -f -15 > "$viruses_folder"/final_outputs/putative_vir_seq_metrics_summary/quality_checkv_summary.tsv

#Copy Taxa estimation file summary
mv -f $viruses_folder'/taxa_estimation_vir_vcontact2_summary.csv' "$viruses_folder"/final_outputs/putative_vir_seq_metrics_summary/taxa_estimation_vir_vcontact2_summary.csv

#Gene
yes | cp "$viruses_folder"/investigation/vibrant/VIBRANT_final_assembly/VIBRANT_results_final_assembly/VIBRANT_genbank_table_final_assembly.tsv "$viruses_folder"/final_outputs/putative_vir_seq_metrics_summary/vir_contigs_genbank_annotation_vibrant.tsv

yes | cp "$viruses_folder"/investigation/vibrant/VIBRANT_final_assembly/VIBRANT_results_final_assembly/VIBRANT_summary_results_final_assembly.tsv "$viruses_folder"/final_outputs/putative_vir_seq_metrics_summary/vir_contigs_annotation_vibrant_summary.tsv
