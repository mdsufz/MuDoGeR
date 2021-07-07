# TO DO STILL:
# FIX ENVIRONMENT WITH VCONTACT INSTALLED INSIDE
# ADJUST PATHS OF TOOLS CALLED IN THE SCRIPT BELOW

conda activate /gpfs1/data/msb/tools/vcontact2/Vontact2

uvigs="$1"
folder_taxonomy="$2"
cores="$3"

date
# creates output directory
mkdir $2

# runs prodigal to generate 1st input to vcontact
prodigal -i "$uvigs" -o  "$folder_taxonomy"/viral_genomes.genes -a "$folder_taxonomy"/viral_genomes.faa -p meta

# runs auxiliary script to generate 2nd input to vcontact
vcontact2_gene2genome -p "$folder_taxonomy"/viral_genomes.faa -o "$folder_taxonomy"/viral_genomes_g2g.csv -s 'Prodigal-FAA'

# runs vcontact
vcontact -t "$cores" --raw-proteins "$folder_taxonomy"/viral_genomes.faa --rel-mode 'Diamond' --proteins-fp "$folder_taxonomy"/viral_genomes_g2g.csv --db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gpfs1/data/msb/tools/vcontact2/cluster_one-1.0.jar --output-dir "$2"/vcontact-output


# process output of vcontact (TO BE DEFINED)
#cat "$folder_taxonomy"/vcontact-output/genome_by_genome_overview.csv  | grep NODE    | cut -f1 -d',' > "$folder_taxonomy"/AUX-1
#while read l; do   grep "$l" "$2"/genome_by_genome_overview.csv | tr '\' ',' | cut -f1,2,19 -d','  ; done < "$2"/AUX-1 > "$2"/AUX-2
