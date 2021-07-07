
#module load GCCcore/8.3.0 Python/3.7.4 Anaconda3/5.3.0
conda activate /gpfs1/data/msb/tools/vcontact2/Vontact2
module load Java/11.0.2



################## COMMAND 

# $1 is the viral genomes
# $2 is the outputdir


date
# creates output directory
mkdir $2

# runs prodigal to generate 1st input to vcontact
/home/brizolat/rodolfo/softwares/dastool/Prodigal/prodigal -i $1 -o  "$2"/viral_genomes.genes -a "$2"/viral_genomes.faa -p meta

# runs auxiliary script to generate 2nd input to vcontact
/gpfs1/data/msb/tools/vcontact2/Vontact2/bin/vcontact2_gene2genome -p "$2"/viral_genomes.faa -o "$2"/viral_genomes_g2g.csv -s 'Prodigal-FAA'

# runs vcontact
/gpfs1/data/msb/tools/vcontact2/Vontact2/bin/vcontact -t ${NSLOTS:-1} --raw-proteins "$2"/viral_genomes.faa --rel-mode 'Diamond' --proteins-fp "$2"/viral_genomes_g2g.csv --db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gpfs1/data/msb/tools/vcontact2/cluster_one-1.0.jar --output-dir "$2"/vcontact-output


# process output of vcontact
cat "$2"/vcontact-output/genome_by_genome_overview.csv  | grep NODE    | cut -f1 -d',' > "$2"/AUX-1

while read l; do   grep "$l" "$2"/genome_by_genome_overview.csv | tr '\n' ',' | cut -f1,2,19 -d','  ; done < "$2"/AUX-1 > "$2"/AUX-2
