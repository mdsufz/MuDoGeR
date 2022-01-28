#!/bin/bash

# TO DO STILL:
# FIX ENVIRONMENT WITH VCONTACT INSTALLED INSIDE
# ADJUST PATHS OF TOOLS CALLED IN THE SCRIPT BELOW


conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database


uvigs="$1"
folder_taxonomy="$2"
cores="$3"

date
# creates output directory
mkdir -p $2

conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/vcontact2_env
# runs prodigal to generate 1st input to vcontact
if [ ! -f "$folder_taxonomy"/viral_genomes.genes ];
then
echo "-> step 1/4"
prodigal -i "$uvigs" -o  "$folder_taxonomy"/viral_genomes.genes -a "$folder_taxonomy"/viral_genomes.faa -p meta
else :
fi


# runs auxiliary script to generate 2nd input to vcontact
if [ ! -f "$folder_taxonomy"/viral_genomes_g2g.csv ];
then
echo "-> step 2/4"
vcontact2_gene2genome -p "$folder_taxonomy"/viral_genomes.faa -o "$folder_taxonomy"/viral_genomes_g2g.csv -s 'Prodigal-FAA'
else :
fi

# runs vcontact
if [ !  -f "$folder_taxonomy"/vcontact-output/genome_by_genome_overview.csv ];
then
echo "-> step 3/4"
cluster_1_path="$(echo "$(echo $PATH | cut -f1 -d':')"/cluster_one-1.0.jar)"
vcontact2 -t "$cores" --raw-proteins "$folder_taxonomy"/viral_genomes.faa --rel-mode 'Diamond' --proteins-fp "$folder_taxonomy"/viral_genomes_g2g.csv --db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin "$cluster_1_path" --output-dir "$2"/vcontact-output
else :
fi

# process output of vcontact (TO BE DEFINED)
if [ ! -f "$folder_taxonomy"/AUX-2 ];
then
echo "-> step 4/4"
cat "$folder_taxonomy"/vcontact-output/genome_by_genome_overview.csv  | grep length   | cut -f1 -d',' > "$folder_taxonomy"/AUX-1
while read l; do   grep "$l" "$folder_taxonomy"/vcontact-output/genome_by_genome_overview.csv  |  cut -f1,2,19 -d','  ; done < "$folder_taxonomy"/AUX-1 > "$folder_taxonomy"/AUX-2
conda deactivate
else :
fi
