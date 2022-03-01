#!/bin/bash

# this bash script generates some simple genome statistics
# should be used like: N50.sh Multi_fasta_file

##### Base of the script: https://github.com/hcdenbakker/N50.sh/blob/master/N50.sh

output_file="$1"/metrics/prok_genomes_stats.tsv
#echo "genome_name""\t""number_of_scaffolds""\t""largest_scaffold_size""\t""N50"      "\t""N90"

echo -e "genome_name\tgenome_size\tnumber_of_scaffolds\tlargest_scaffold_size\tN50\tN90" > $output_file
for genome in "$1"/binning/unique_bins/*.fa;
do
genome_name="$(echo "$genome" | rev | cut -f1 -d'/' | rev )";
mkdir temp; #check location
#get contig lengths, ordered from large to small; 1. remove newlines in sequences, 2. remove fasta headers
# 3. get contig sizes (line lengths) and order them from large to small.
cat $genome | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'|
sed '/^>/ d'| awk '{ print length($0) }' | sort -gr > temp/contig_lengths.txt;
#number of contigs
Y=$(cat temp/contig_lengths.txt | wc -l);
#sum of contig_lengths
X=$(paste -sd+ temp/contig_lengths.txt | bc);
# cumulative contig lengths
awk 'BEGIN {sum=0} {sum= sum+$0; print sum}' temp/contig_lengths.txt > temp/contig_lengths_cum.txt;
# get cumulative contig contributions (%) to the entire genome
awk -v var=$X 'BEGIN {FS=OFS=","} {for (i=1; i<=NF; i++) $i/=var;}1' temp/contig_lengths_cum.txt > temp/cum_perc.txt; 
# join results
paste temp/contig_lengths.txt temp/cum_perc.txt > temp/matrix.txt;
# get N50, N90, largest scaffold/contig
N50=$(awk '$2 >= 0.50' temp/matrix.txt |head -1| awk '{ print $1}');
N90=$(awk '$2 >= 0.90' temp/matrix.txt |head -1| awk '{ print $1}');
large_contig=$(head -1 temp/contig_lengths.txt);
rm -fr temp;

#Echo "$genome" "\t" "$Y""\t"$X"\t" "$large_contig"       "\t" "$N50" "\t" "$N90"
echo -e "$genome_name\t$X\t$Y\t$large_contig\t$N50\t$N90" >> $output_file
done

### summarize results in one table 

echo -e "OTU\tcompleteness\tcontamination\tstr.heterogeneity\ttaxonomy\tgenome_size\t#scaffolds\tlargest_scaff\tN50\tN90\tprokka_known\tprokka_unknown" > genome_metrics.tsv ;
for d in  binning/unique_bins/*;
do bin="$(echo $d | rev | cut -f1 -d'/' | rev | sed "s/.fa//g")"; 
echo -e "$bin\t\c"; 
tax="$(grep "$bin" metrics/GTDBtk_taxonomy/*.summ* | cut -f2)";
qual="$(grep "$bin" metrics/checkm_qc/outputcheckm.tsv | cut -f12,13,14 )";
echo -e "$qual\t\c";
echo -e "$tax\t\c";
metrics="$(grep "$bin" metrics/prok_genomes_stats.tsv  | cut -f2-10)";
echo -e "$metrics\t\c"; 
prokka_known="$(tail -n +2 metrics/prokka/"$bin"/PROKKA*tsv  | grep -v hypothetical | wc -l )";
prokka_unknown="$(tail -n +2 metrics/prokka/"$bin"/PROKKA*tsv  | grep hypothetical | wc -l )";
echo -e "$prokka_known\t$prokka_unknown";done >> genome_metrics.tsv
