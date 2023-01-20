#!/bin/bash

# this bash script generates some simple genome statistics
# should be used like: N50.sh Multi_fasta_file

##### Base of the script: https://github.com/hcdenbakker/N50.sh/blob/master/N50.sh
# loading conda environment
echo '------- START MODULE 2-8 BIN METRICS'
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

# load bbtools env
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/bbtools_env

mkdir -p "$1"/metrics/genome_statistics

prokaryotes_folder=$1
output_path="$1"/metrics/genome_statistics
output_file="$output_path"/prok_genomes_stats.tsv
cores=$2


#Create file
echo -e "genome_name\tgenome_size\tnumber_of_scaffolds\tlargest_scaffold_size\tN50\tN90" > $output_file

for genome in "$1"/binning/unique_bins/*.fa;
do
genome_name="$(echo "$genome" | rev | cut -f1 -d'/' | rev )";
mkdir -p "$1"/metrics/temp;
temp_path="$1"/metrics/temp;
#get contig lengths, ordered from large to small; 1. remove newlines in sequences, 2. remove fasta headers
# 3. get contig sizes (line lengths) and order them from large to small.
cat $genome | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'|
sed '/^>/ d'| awk '{ print length($0) }' | sort -gr > $temp_path/contig_lengths.txt;
#number of contigs
Y=$(cat $temp_path/contig_lengths.txt | wc -l);
#sum of contig_lengths
X=$(paste -sd+ $temp_path/contig_lengths.txt | bc);
# cumulative contig lengths
awk 'BEGIN {sum=0} {sum= sum+$0; print sum}' $temp_path/contig_lengths.txt > $temp_path/contig_lengths_cum.txt;
# get cumulative contig contributions (%) to the entire genome
awk -v var=$X 'BEGIN {FS=OFS=","} {for (i=1; i<=NF; i++) $i/=var;}1' $temp_path/contig_lengths_cum.txt > $temp_path/cum_perc.txt; 
# join results
paste $temp_path/contig_lengths.txt $temp_path/cum_perc.txt > $temp_path/matrix.txt;
# get N50, N90, largest scaffold/contig
N50=$(awk '$2 >= 0.50' $temp_path/matrix.txt |head -1| awk '{ print $1}');
N90=$(awk '$2 >= 0.90' $temp_path/matrix.txt |head -1| awk '{ print $1}');
large_contig=$(head -1 $temp_path/contig_lengths.txt);


rm -fr $temp_path;

#Echo "$genome" "\t" "$Y""\t"$X"\t" "$large_contig"       "\t" "$N50" "\t" "$N90"
echo -e "$genome_name\t$X\t$Y\t$large_contig\t$N50\t$N90" >> $output_file
done

### Summarize results in one table 

echo -e "OTU\tcompleteness\tcontamination\tstr.heterogeneity\ttaxonomy\tgenome_size\t#scaffolds\tlargest_scaff\tN50\tN90\tprokka_known\tprokka_unknown" > $output_path/genome_metrics.tsv ;
for d in  $1/binning/unique_bins/*;
do bin="$(echo $d | rev | cut -f1 -d'/' | rev | sed "s/.fa//g")"; 
echo -e "$bin\t\c"; 
tax="$(grep "$bin" $1/metrics/GTDBtk_taxonomy/*.summ* | cut -f2)";
qual="$(grep "$bin" $1/metrics/checkm_qc/outputcheckm.tsv | cut -f12,13,14 )";
echo -e "$qual\t\c";
echo -e "$tax\t\c";
metrics="$(grep "$bin" $output_path/prok_genomes_stats.tsv  | cut -f2-10)";
echo -e "$metrics\t\c"; 
prokka_known="$(tail -n +2 $1/metrics/prokka/"$bin"/PROKKA*tsv  | grep -v hypothetical | wc -l )";
prokka_unknown="$(tail -n +2 $1/metrics/prokka/"$bin"/PROKKA*tsv  | grep hypothetical | wc -l )";
echo -e "$prokka_known\t$prokka_unknown"; done >> $output_path/genome_metrics.tsv


##Run BBTOOLS on selected MAGs
cd "$1"/binning/unique_bins/

"$MUDOGER_DEPENDENCIES_ENVS_PATH"/bbtools_env/bin/statswrapper.sh *.fa > $output_path/bbtools.tsv

## Filter Good quality bins (Complet - 5*Contamination >= 50) based on CheckM
cat $output_path/genome_metrics.tsv | awk 'BEGIN {FS="\t"};NR==1; ($2 - (5*$3)) >= 50' > "$1"/MAGS_results.tsv

cd -

## Colect all relevant results in one final output folder

mkdir -p "$prokaryotes_folder"/final_outputs
mkdir -p "$prokaryotes_folder"/final_outputs/all_bins_seq
mkdir -p "$prokaryotes_folder"/final_outputs/only_mags_seq
mkdir -p "$prokaryotes_folder"/final_outputs/bins_metrics_summary
mkdir -p "$prokaryotes_folder"/final_outputs/bins_genes_prokka_summary

#copy bins
yes | cp "$prokaryotes_folder"/binning/unique_bins/*.fa "$prokaryotes_folder"/final_outputs/all_bins_seq/
#Copy taxa
yes | cp "$prokaryotes_folder"/metrics/GTDBtk_taxonomy/*bac*.summ* "$prokaryotes_folder"/final_outputs/bins_metrics_summary/taxa_bins_gtdbtk_summary.tsv
yes | cp "$prokaryotes_folder"/metrics/GTDBtk_taxonomy/*ar*.summ* "$prokaryotes_folder"/final_outputs/bins_metrics_summary/arc_taxa_bins_gtdbtk_summary.tsv
#Copy quality
yes | cp "$prokaryotes_folder"/metrics/checkm_qc/outputcheckm.tsv "$prokaryotes_folder"/final_outputs/bins_metrics_summary/qual_bins_checkm_summary.tsv
#copy gene annotation
for d in  $1/binning/unique_bins/*;
do bin="$(echo $d | rev | cut -f1 -d'/' | rev | sed "s/.fa//g")";
yes | cp "$prokaryotes_folder"/metrics/prokka/"$bin"/PROKKA*tsv "$prokaryotes_folder"/final_outputs/bins_genes_prokka_summary/"$bin"_genes_prokka.tsv
done

#Copy complete summary
yes | cp $output_path/genome_metrics.tsv "$prokaryotes_folder"/final_outputs/allbins_metrics_summary.tsv

#Move mags summary
mv "$prokaryotes_folder"/MAGS_results.tsv "$prokaryotes_folder"/final_outputs/mags_results_summary.tsv

#Copy only MAGs seq
cat "$prokaryotes_folder"/final_outputs/mags_results_summary.tsv | cut -f1 | tail -n+2 > "$prokaryotes_folder"/final_outputs/tmp
while read mag; do yes | cp "$prokaryotes_folder"/binning/unique_bins/$mag.fa "$prokaryotes_folder"/final_outputs/only_mags_seq;
done < "$prokaryotes_folder"/final_outputs/tmp

rm -f "$prokaryotes_folder"/final_outputs/tmp

conda deactivate




