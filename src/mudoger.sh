
kmer



######$ -N khmertest
#$ -l h_rt=3:00:00
#$ -l h_vmem=10G
#$ -binding linear:1



 # output files
#$ -o /gpfs1/data/msb/thym/mudoger_tutorial/scripts/LOGS/$JOB_NAME-$JOB_ID.out
#$ -e /gpfs1/data/msb/thym/mudoger_tutorial/scripts/LOGS/$JOB_NAME-$JOB_ID.err



input_file=$1
output_file33=$2
output_file55=$3
k_size33=33
k_size55=55


echo "memory is" "$memory"
echo "started khmer"
date

#khmer count
/data/msb/tools/miniconda/miniconda2/bin/unique-kmers.py -k $k_size33 -R $output_file33 $input_file

/data/msb/tools/miniconda/miniconda2/bin/unique-kmers.py -k $k_size55 -R $output_file55 $input_file


echo "end counting"
date


echo -e "library_id\t\c" > input.tsv;   echo -e "$(head -n1 $output_file33  | cut -f1 -d' ')\t\c" >> input.tsv;  echo  "$(head -n1 $output_file55 | cut -f1 -d' ')" >> input.tsv

#Load modules
module load   GCC/8.3.0  OpenMPI/3.1.4

module load R/3.6.2-2


Rscript predict.R input.tsv

cat metaspades_prediction.tsv  | rev | cut -f1 | rev | grep -v max
sed 's/\t/,/g'  metaspades_prediction.tsv > metaspades_prediction.csv

### PROK
# resoruces
#$ -l h_rt=40:00:00
#$ -l h_vmem=60G
###$ -binding linear:1  # for single core       #
#$ -pe smp 4-28         # for multi-core




# output files
#$ -o /gpfs1/data/msb/thym/mudoger_tutorial/scripts/LOGS/$JOB_NAME-$JOB_ID.out
#$ -e /gpfs1/data/msb/thym/mudoger_tutorial/scripts/LOGS/$JOB_NAME-$JOB_ID.err



# loading modules
module load  GCC/6.4.0-2.28  OpenMPI/2.1.2 Python/2.7.14


source /data/msb/tools/miniconda/miniconda2/bin/activate
conda activate env_metawrap_v1.2.3

echo "startd binning"
date
mkdir $1
prok="$1"

bin="$prok"


# metawrap READ_QC
metawrap  binning -o $bin -t ${NSLOTS:-1} -a $2 --run-checkm --metabat2 --maxbin2 --concoct  $3 $4
#metawrap  binning -o $1 -t ${NSLOTS:-1} -a $2  --metabat2 --maxbin2 --concoct  $3 $4

touch "$1"/"finished_binning"

echo "end binning"
date


# loading modules
module load gcc/4/8.1-3

   source /data/msb/tools/miniconda/miniconda2/bin/activate
#conda activate /data/msb/tools/metawrap/metawrap_env
#conda activate /data/msb/tools/metawrap/THIRD_INSTALL/metawrap_env_3
conda activate env_metawrap_v1.2.3

echo "startd bin refinement"
date


mkdir $prok/bact_ref
outb="$prok/bact_ref"
mkdir $prok/arch_ref
outa="$prok/arch_ref"

con="$prok/concoct_bins"
met="$prok/metabat2_bins"
max="$prok/maxbin2_bins"
# metawrap READ_QC

# metawrap read_qc --skip-bmtagger -1 A_5_V_60-Empirical_fw.fastq -2 A_5_V_60-Empirical_bw.fastq -t 2 -o READ_QC/A5V60_QCed.fq
#1
metawrap bin_refinement -o "$outb" -t "${NSLOTS:-1}" -A "$con" -B "$met" -C "$max" -c 50 -x 10


metawrap bin_refinement -o "$outa" -t "${NSLOTS:-1}" -A "$con" -B "$met" -C "$max" -c 40 -x 30


touch "$1"/refinement_finished

echo "end bin refinement"
date

ch_inpb="$outb/metawrap_50_10_bins"
ch_inpa="$outa/metawrap_40_30_bins"





# loading modules

module load   GCC/6.4.0-2.28  OpenMPI/2.1.2 Python/3.6.4

  source /gpfs1/data/msb/tools/GTDB/gtdbtk-v1.3.0/bin/activate

# necessary path variable
export GTDBTK_DATA_PATH=/gpfs1/data/msb/tools/GTDB/external_data/release95

mkdir $prok/bacteria_output_tax_directory
mkdir $prok/archaea_output_tax_directory





gtdbtk  classify_wf --extension  fa  --cpus ${NSLOTS:-1} --genome_dir $ch_inpb  --out_dir $prok/bacteria_output_tax_directory

gtdbtk  classify_wf --extension  fa  --cpus ${NSLOTS:-1} --genome_dir $ch_inpa  --out_dir $prok/archaea_output_tax_directory


cd $prok/bacteria_output_tax_directory/classify/intermediate_results
sed 's/\t/,/g' gtdbtk.bac120.classification_pplacer.tsv > gtdbtk.bac120.classification_pplacer.csv
cd ..
cd ..
cd ..
cd ..
cd $prok/archaea_output_tax_directory/classify/intermediate_results
sed 's/\t/,/g' gtdbtk.bac120.classification_pplacer.tsv > gtdbtk.bac120.classification_pplacer.csv

cd ..
cd ..
cd ..
cd ..


# loading modules

module load GCCcore/8.3.0 Python/3.7.4


task=""

date
echo "starting "$task

source /data/msb/tools/checkm/new_install/checkm/bin/activate


#command
checkm lineage_wf -t ${NSLOTS:-1}  --reduced_tree --tab_table -x fa -f $prok/checkm_archaea/checkm.tsv  $ch_inpa $prok/checkm_archaea

#command
checkm lineage_wf -t ${NSLOTS:-1}  --reduced_tree --tab_table -x fa -f $prok/checkm_bacteria/checkm.tsv $ch_inpb $prok/checkm_bacteria


date
echo "finishing "$task

#FILTERING CHECKM results

if [ ! -z $5 ]
then
        q="$5"
        cat $prok/checkm_archaea/checkm.tsv | cut -f1,12,13 | awk -F '\t' '{ if ( $2- 5* $3 >= '$q' ) { print } }' > $prok/checkm_archaea/filtered_checkm.tsv

        cat $prok/checkm_bacteria/checkm.tsv | cut -f1,12,13 | awk -F '\t' '{ if ( $2- 5* $3 >= '$q' ) { print } }' > $prok/checkm_bacteria/filtered_checkm.tsv

else
        echo "no extra filtering is chosen"
fi


cd $prok/checkm_archaea

sed 's/\t/,/g' checkm.tsv > checkm.csv
sed 's/\t/,/g' filtered_checkm.tsv > filtered_checkm.csv

cd ..
cd ..

cd  $prok/checkm_bacteria
sed 's/\t/,/g' checkm.tsv >  checkm.csv
sed 's/\t/,/g' filtered_checkm.tsv > filtered_checkm.csv
cd ..
cd ..


module purge
conda deactivate
# loading modules
module load gcc/4/8.1-3

source /data/msb/tools/miniconda/miniconda2/bin/activate
conda activate /data/msb/tools/metawrap/metawrap_env

unset LD_PRELOAD


mkdir $prok/Annotation_archaea
mkdir $prok/Annotation_bacteria

metawrap annotate_bins -o $prok/Annotation_archaea -t "${NSLOTS:-1}" -b $ch_inpa

metawrap annotate_bins -o $prok/Annotation_bacteria -t "${NSLOTS:-1}" -b $ch_inpb


#### VIRAL MODULE


#$ -l h_rt=50:00:00
#$ -l h_vmem=30G
####$ -binding linear:1         # for single core       #
#$ -pe smp 2-28         # for multi-core



# output files
#$ -o /gpfs1/data/msb/thym/mudoger_tutorial/scripts/LOGS/$JOB_NAME-$JOB_ID.out
#$ -e /gpfs1/data/msb/thym/mudoger_tutorial/scripts/LOGS/$JOB_NAME-$JOB_ID.err



# STEP 1 RECOVERY USING VIBRANT
# modules to load
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 Python/3.6.6
# environment to load
source /data/msb/tools/vibrant/env_vibrant_v1.2.1/bin/activate

mkdir $1
output_viral="$1"

mkdir $output_viral/vibrant_folder
cd $output_viral/vibrant_folder


python /data/msb/tools/vibrant/VIBRANT/VIBRANT_run.py -i $2


module purge
conda deactivate

cd ..
cd ..


module purge
conda deactivate


  #STEP 2 VIRFINDER RECOVERY

# loading modules 
#module load  GCC/7.3.0-2.30  OpenMPI/3.1.1 R/3.5.1 
module load  GCC/6.4.0-2.28  OpenMPI/2.1.2 R/3.4.4-X11-20180131

echo "started virfinder"
date

#mkdir $1
#output_viral="$1"
mkdir $output_viral/virfinder_folder
inp_virfinder=$output_viral/virfinder_folder
#command
Rscript /data/msb/tools/virfinder/virfinder_script.r $inp_virfinder $2 $inp_virfinder/virfinder.tsv

echo "end virfinder"
date


# STEP 3 VIRSORTER RECOVERY 


# loading modules
module load gcc/4/8.1-3


#conda environments
source /data/msb/tools/miniconda/miniconda2/bin/activate
conda activate /data/msb/tools/virsorter/virsorter_env



echo "started virsorter"
virs_inp=$output_viral/virsorter_folder

#command
wrapper_phage_contigs_sorter_iPlant.pl -f $2 --wdir $virs_inp


echo "end virsorter"
module purge
conda deactivate
                      


#STEP 4 FILTERING OF THE OUTPUTS AND EXTRACTION OF SEQUENCES 

#Formation of destination file for the filtered data
mkdir $output_viral/viral_filtered_data


# STEP 4.1 FILTERING VIRFINDER.
#Looking for sequences with q-value <= 0.01 and length >= 1000 


#Command 
cat $inp_virfinder/virfinder.tsv | awk -F'\t' '{ if ( $4 <= 0.01) print }' | awk -F'_' '{ if ( $4 >= 1000) print  }' | cut -f2 | sed "s/\"//g" > $output_viral/viral_filtered_data/virfinder_filtered_data


#STEP 4.2 FILTERING VIRSORTER 

#Looking only for sequences from categories 1 and 2

#Command
cat $virs_inp/Predicted_viral_sequences/VIRSorter_cat-{1..2}*fasta | grep ">" | sed "s/>VIRSorter_//g"  | sed "s/-cat_2//g" | sed "s/-cat_1//g" | sed 's/\(.*\)_/\1./' > $output_viral/viral_filtered_data/virsorter_filtered_data

#STEP 4.3 FILTERING VIBRANT 

vibrant_filt_inp="$output_viral/vibrant_folder/VIBRANT_P19013_NLE6_assembly/VIBRANT_phages_P19013_NLE6_assembly"

#Command 
cat $vibrant_filt_inp/P19013_NLE6_assembly.phages_combined.fna | grep ">" | sed "s/_fragment_1//g;s/>//g" > $output_viral/viral_filtered_data/VIBRANT_filtere_data


#STEP 4.4 COMBINATION OF THE FILTERED OUTPUTS AND REMOVAL OF THE REPEATED SEQUENCES


cd $output_viral/viral_filtered_data
#Command

cat * virfinder_filtered_data virsorter_filtered_data VIBRANT_filtere_data | sort | uniq > COMBINED_VIRAL_PARTICLES_FOR_EXTRACTION

cd ..
cd ..


mkdir $output_viral/extracted_and_dereplication_fasta

conda deactivate
module purge

#command

python /data/msb/thym/mudoger_tutorial/viral_test_folder/extract_fa.py $output_viral/viral_filtered_data/COMBINED_VIRAL_PARTICLES_FOR_EXTRACTION $2 $output_viral/extracted_and_dereplication_fasta/VIRAL_PARTICLES.fa

#STEP 5 DEREPLICATION 

derep_input="$output_viral/extracted_and_dereplication_fasta/VIRAL_PARTICLES.fa"


#Removal of new modules
module purge

task=""
date
echo "starting "$task


#Command 
/gpfs1/data/msb/tools/MUMmer/new_install/stampede-clustergenomes/bin/Cluster_genomes.pl -f $derep_input  -c 70 -i 95


date
echo "finishing "$task



#Cleaning the output directory 
echo "Cleaning the directory"
cd $output_viral
rm -rf viral_filtered_data
mv  extracted_and_dereplication_fasta dereplication_folder
cd dereplication_folder
rm -rf VIRAL_PARTICLES.fa
cd ..
mv $output_viral/dereplication_folder/VIRAL_PARTICLES_95-70.fna $output_viral/dereplication_folder/VIRAL_PARTICLES_95-70.fa



mkdir $output_viral/initial_recovery_folder
mv $output_viral/vibrant_folder $output_viral/initial_recovery_folder
mv $output_viral/virsorter_folder $output_viral/initial_recovery_folder
mv $output_viral/virfinder_folder $output_viral/initial_recovery_folder


module purge


#CHECKv

# loading modules
module Anaconda3/5.3.0


conda activate /data/msb/tools/vcheck/conda-vcheck-env
export CHECKVDB=/data/msb/tools/vcheck/checkv-db-v0.6


# INPUTS
#1 is a fasta file containing all viral particles (after dereplication with stampede clustergenomes) 
#2 is the output directory where the results will be dumped (this file will be created)
# USAGE

new_der="$output_viral/dereplication_folder/VIRAL_PARTICLES_95-70.fa"




# start
checkv end_to_end  $new_der $output_viral/quality_folder  -t ${NSLOTS:-1}

date
echo "finishing "$task

module purge
conda deactivate

# loading modules
module load GCCcore/8.3.0 Python/3.7.4 Anaconda3/5.3.0

conda activate /gpfs1/data/msb/tools/vcontact2/Vontact2

module load Java/11.0.2



# INPUTS
#1 is a fasta file containing all viral particles (after dereplication with stampede clustergenomes) 
#2 is the output directory where the results will be dumped (this file will be created)
# USAGE
### qsub -m beas -M <my ufz e-mail> -N <jobname> sub_script.sh  <input_1> <input_2> 



# $1 is the viral genomes to be investigated. One fasta file is provided.
# $2 is the outputdir. It will be created.
#start
# creates output directory
mkdir  $output_viral/taxonomy_folder
tax_out=$output_viral/taxonomy_folder
# runs prodigal to generate 1st input to vcontact
/gpfs1/data/msb/tools/mudoger/miniconda2/bin/prodigal -i $new_der -o  "$tax_out"/viral_genomes.genes -a "$tax_out"/viral_genomes.faa -p meta

# runs auxiliary script to generate 2nd input to vcontact
vcontact2_gene2genome -p "$tax_out"/viral_genomes.faa -o "$tax_out"/viral_genomes_g2g.csv -s 'Prodigal-FAA'

# runs vcontact
vcontact -t ${NSLOTS:-1} --raw-proteins "$tax_out"/viral_genomes.faa --rel-mode 'Diamond' --proteins-fp "$tax_out"/viral_genomes_g2g.csv --db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gpfs1/data/msb/tools/vcontact2/cluster_one-1.0.jar --output-dir "$tax_out"/vcontact-output
 
 # process output of vcontact
cat "$tax_out"/vcontact-output/genome_by_genome_overview.csv  | grep NODE    | cut -f1 -d',' > "$tax_out"/AUX-1

while read l; do   grep "$l", genome_by_genome_overview.csv | tr '\n' ',' | cut -f1,2,19 -d','  ; done < "$tax_out"/AUX-1 > "$tax_out"/AUX-2

## This is for big datasets

#mkdir  $output_viral/taxonomy_folder
#tax_out=$output_viral/taxonomy_folder
# runs prodigal to generate 1st input to vcontact
#/gpfs1/data/msb/tools/mudoger/miniconda2/bin/prodigal -i $output_viral/quality_folder/cleaned_contigs.fna -o "$tax_out"/viral_genomes.genes -a "$tax_out"/viral_genomes.faa -p meta

# runs auxiliary script to generate 2nd input to vcontact
#vcontact2_gene2genome -p "$tax_out"/viral_genomes.faa -o "$tax_out"/viral_genomes_g2g.csv -s 'Prodigal-FAA'

# runs vcontact
#vcontact -t ${NSLOTS:-1} --raw-proteins "$tax_out"/viral_genomes.faa --rel-mode 'Diamond' --proteins-fp "$tax_out"/viral_genomes_g2g.csv --db 'ProkaryoticViralRefSeq94-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /gpfs1/data/msb/tools/vcontact2/cluster_one-1.0.jar --output-dir "$tax_out"/vcontact-output


# process output of vcontact
#cat "$tax_out"/vcontact-output/genome_by_genome_overview.csv  | grep NODE    | cut -f1 -d',' > "$tax_out"/AUX-1

#while read l; do   grep "$l", genome_by_genome_overview.csv | tr '\n' ',' | cut -f1,2,19 -d','  ; done < "$tax_out"/AUX-1 > "$tax_out"/AUX-2



module purge
conda deactivate


 ## WISH

# 

# INPUTS
#1 is a fasta file containing all viral particles (after dereplication with stampede clustergenomes) 
#2 is the directory containg all prokaryotic bins (the potential host genomes)
#3 is the output directory where the results will be dumped (this file will be created)




if [ ! -z $3 ]
then

        mkdir "$1/wish_folder"
        wish="$1/wish_folder"
        mkdir "$wish"/potential_host_genomes
        mkdir "$wish"/viral_particles
        mkdir "$wish"/output_results
        #2 cp and prepare data
        cp "$3"/*fa "$wish"/potential_host_genomes
        python /data/msb/tools/wish/split-all-seq.py "$new_der" "$wish"/viral_particles/viral-particle
        echo "steps: 2/5"
        #3 build model
        /data/msb/tools/wish/WIsH/WIsH -c build -g "$wish"/potential_host_genomes/ -m "$wish"/modelDir
        echo "steps: 3/5"
        #4 run prediction
        /data/msb/tools/wish/WIsH/WIsH -c predict -g "$wish"/viral_particles/ -m "$wish"/modelDir/ -r "$wish"/output_results/ -b
        echo "steps: 4/5"
        #5 convert to p-value
        /data/msb/tools/wish/WIsH/WIsH -c predict -g "$wish"/viral_particles/ -m "$wish"/modelDir/ -r  "$wish"/output_results/ -b -p -n /data/msb/tools/wish/WIsH/KeggGaussianFits.tsv
        echo "steps: 5/5"
else
        echo "No host identification"
fi


# THis is in case of large datasets

#if [ ! -z $3 ]
#then

 #       mkdir "$1/wish_folder"
  #      wish="$1/wish_folder"
   #     mkdir "$wish"/potential_host_genomes
    #    mkdir "$wish"/viral_particles
     #   mkdir "$wish"/output_results
        #2 cp and prepare data
      #  cp "$3"/*fa "$wish"/potential_host_genomes
       # python /data/msb/tools/wish/split-all-seq.py "$output_viral/quality_folder/cleaned_contigs.fna" "$wish"/viral_particles/viral-particle
       # echo "steps: 2/5"
        #3 build model
       # /data/msb/tools/wish/WIsH/WIsH -c build -g "$wish"/potential_host_genomes/ -m "$wish"/modelDir
       # echo "steps: 3/5"
        #4 run prediction
       # /data/msb/tools/wish/WIsH/WIsH -c predict -g "$wish"/viral_particles/ -m "$wish"/modelDir/ -r "$wish"/output_results/ -b
       # echo "steps: 4/5"
        #5 convert to p-value
       # /data/msb/tools/wish/WIsH/WIsH -c predict -g "$wish"/viral_particles/ -m "$wish"/modelDir/ -r  "$wish"/output_results/ -b -p -n /data/msb/tools/wish/WIsH/KeggGaussianFits.tsv
       # echo "steps: 5/5"
#else
 #       echo "No host identification"
#fi







cd $output_viral/initial_recovey_folder/virfinder_folder
sed 's/\t/,/g' virfinder.tsv > virfinder.csv
cd ..
cd ..
cd ..

cd $output_viral/quality_folder
sed 's/\t/,/g' completeness.tsv > completeness.csv
sed 's/\t/,/g' contamination.tsv > contamination.csv
sed 's/\t/,/g' quality_summary.tsv > quality_summary.csv
                                                           
cd ..
cd ..




echo "VIRAL RECOVERY IS FINISHED!"

### euk module


#echo "First arg: $1"
#echo "Second arg: $2"
#echo "Third arg: $3"
#echo "Fourth arg: $4"


module load GCC/6.4.0-2.28 Python/3.6.4


mkdir $2
mkdir $3


source /data/msb/tools/eukrep/env_p3.6.4/bin/activate


echo "beginning eukrep"
date

EukRep -i $1 --prokarya $2/prok.fasta -o $3/euk.fasta


echo "finishing eukrep"
date




# loading modules
module load gcc/4/8.1-3


source /data/msb/tools/miniconda/miniconda2/bin/activate
conda activate /data/msb/tools/metawrap/metawrap_env


#export PERL5LIB=/data/msb/tools/metawrap/metawrap_env/lib/perl5/site_perl/5.22.2/

echo "startd binning"
date

mkdir $3/conc
con="$3/conc"

conbin="$3/euk.fasta"

# metawrap READ_QC
metawrap  binning -o $con -t ${NSLOTS:-1} -a $conbin --concoct  $4 $5

echo "end binning"
mkdir $con/euk_concoct_bins

for d in "$con"/concoct_bins/*; do id=$(echo $d | cut -f8 -d'/'); nb=$con/euk_concoct_bins/"$(echo $d | sed "s/bin\./"$id"-ebin./g" | rev | cut -f1 -d'/' | rev )"; echo copying "$d" "$nb" ; cp "$d" "$nb" ;done

cd $con/euk_concoct_bins
find . -size -2500k -delete

cd ..
cd ..
cd ..


echo "End of the first part of eukaryotic module"

#loading modules
module load foss/2019b Perl/5.30.0


#$ -l h_rt=10:00:00
#$ -l h_vmem=10G
#$ -binding linear:1    # for single core       #
###$ -pe smp 15-32      # for multi-core

# output files
#$ -o /data/msb/thym/mudoger_tutorial/scripts/LOGS/$JOB_NAME-$JOB_ID.out
#$ -e /data/msb/thym/mudoger_tutorial/scripts/LOGS/$JOB_NAME-$JOB_ID.err

task="genemark-ES"

date
echo "starting "$task

mkdir "$2"

cp "$1" "$2"

cd "$2"

bin="$(echo $1 | rev | cut -f1 -d'/' | rev )"

perl -Mlocal::lib=/gpfs1/data/msb/tools/genemark/NEW_INSTALL/perl_mods /gpfs1/data/msb/tools/genemark/NEW_INSTALL/gmes_linux_64/gmes_petap.pl  --ES -min_contig 3000 --sequence $2/"$bin"


#MAKER

#loading modules

module load MAKER/2.31.10-foss-2019b-1




echo "maker"
date

genemark_gmhmm_file="$2/output/gmhmm.mod"
genome_fasta="$1"
outputfolder="$2/maker"


mkdir $outputfolder
cp /data/msb/tools/maker2/installation_tom_module/DEFAULT_CONTROL_FILES/* $outputfolder
cp $genome_fasta $outputfolder

cd $outputfolder

sed -i  's|gmhmm=|gmhmm='$genemark_gmhmm_file'|g'  "$outputfolder"/maker_opts.ctl

maker -g $genome_fasta -c ${NSLOTS:-1}

cd *output

fasta_merge -d *_master_datastore_index.log -o OUTPUT

echo "end maker"



#modules BUSCO
 module load GCC/6.4.0-2.28 OpenMPI/2.1.2 Python/3.6.4


source /gpfs0/global/apps/msb/busco/busco/busco_env/bin/activate

echo "start busco"

date

mkdir $outputfolder/busco
bus=$outputfolder/busco
output="$(echo $bus | rev | cut -f1 -d'/' | rev )"

/global/apps/msb/busco/busco/scripts/run_BUSCO.py -i $1 -l /gpfs0/global/apps/msb/busco/eukaryota_odb9 -m prot -o $output

cd $bus

sed 's/'$'\t''/,/g' full_table_fbusco.tsv > full_table_fbusco.csv
cd ..
cd ..



 echo "end busco"

# eukcc
# loading modules
module load Anaconda3/2019.10
module load GCCcore/8.3.0 Python/3.7.4
module load foss/2019b Perl/5.30.0


k="eukcc"

date
echo "starting "$task


source activate /data/msb/tools/eukcc/eukcc_env


# eukcc starts after filtering

#1 is the input eMAG
#2 is the output folder

mkdir $2/eukcc

#command
eukcc --db /gpfs1/data/msb/tools/eukcc/eukccdb  --ncores "${NSLOTS:-1}"   --ncorespplacer 1  --outdir "$2/eukcc" "$1"





date
echo "finishing "$task


               
 








 

                                                                                                                                                                                 




                                                                                                                                                 2,0-1         


                                                                                                                                                     1,1           


