#!/bin/bash

# this scripts asks for the desired location of installation for the databases.
# given the user input, the config file will be edited and all databases via wget, curl, etc and zipping

database_location="$1"

mkdir -p "$database_location"
conda activate mudoger_env
config_file="$(which config.sh)"
source "$config_file"
echo DATABASES_LOCATION="$database_location" > ${config_file/config/database}


############################################### PROKARYOTES ###############################################
### CheckM
mkdir -p "$database_location"/checkm
cd  "$database_location"/checkm
if [ ! -f selected_marker_sets.tsv ]; then
echo 'installing checkm database ...'
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf checkm_data_2015_01_16.tar.gz
rm -fr checkm_data_2015_01_16.tar.gz

CHECKM_DB="$database_location"/checkm #Fixed? we need to test
echo CHECKM_DB="$CHECKM_DB" >> "$config_path"
else echo "-> your CheckM database is ready"
fi


### GTDB-tk
mkdir -p  "$database_location"/"gtdbtk"
cd "$database_location"/"gtdbtk"
if [ ! -d release*  ]; then
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
rm -fr gtdbtk_data.tar.gz

else echo "-> your GTDBtk database is ready"
fi



############################################### VIRUSES ###############################################
## VIBRANT
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
pip install pickle-mixin --quiet
VIBRANT_DB_DIR=$database_location/vibrant
if [ ! -f $VIBRANT_DB_DIR/Pfam-A_v32.HMM.h3p ] ;
then
#echo 'let us download '$VIBRANT_DB_DIR
if [  ! -f $VIBRANT_DB_DIR/vog.hmm.tar.gz ] ; then
wget http://fileshare.csb.univie.ac.at/vog/vog94/vog.hmm.tar.gz -P $VIBRANT_DB_DIR
else :; fi
if [ ! -f $VIBRANT_DB_DIR/Pfam-A.hmm.gz ]; then
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz -P $VIBRANT_DB_DIR
else :; fi
if [ ! -f $VIBRANT_DB_DIR/profiles.tar.gz ]; then
wget ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/profiles.tar.gz -P $VIBRANT_DB_DIR
else :; fi
#echo 'downloaded'
tar -xzf $VIBRANT_DB_DIR/vog.hmm.tar.gz -C $VIBRANT_DB_DIR
gunzip $VIBRANT_DB_DIR/Pfam-A.hmm.gz -d $VIBRANT_DB_DIR
tar -xzf $VIBRANT_DB_DIR/profiles.tar.gz -C $VIBRANT_DB_DIR
#echo 'decompressed'
for v in $VIBRANT_DB_DIR/VOG*.hmm;
do cat $v >> $VIBRANT_DB_DIR/vog_temp.HMM;
done
for k in $VIBRANT_DB_DIR/profiles/K*.hmm;
do cat $k >> $VIBRANT_DB_DIR/kegg_temp.HMM;
done
#echo 'concatenated'
rm -f $VIBRANT_DB_DIR/VOG0*.hmm
rm -f  $VIBRANT_DB_DIR/VOG1*.hmm
rm -f  $VIBRANT_DB_DIR/VOG2*.hmm
rm -rf $VIBRANT_DB_DIR/profiles
#echo 'clean'
prof_names="$(echo $PATH | sed "s/:/\n/g" | grep vibrant | sed "s/bin/share\/vibrant-1.2.0\/databases\/profile_names/g")"
cp -r $prof_names $VIBRANT_DB_DIR
hmmfetch -o $VIBRANT_DB_DIR/VOGDB94_phage.HMM -f $VIBRANT_DB_DIR/vog_temp.HMM $VIBRANT_DB_DIR/profile_names/VIBRANT_vog_profiles.txt
hmmfetch -o $VIBRANT_DB_DIR/KEGG_profiles_prokaryotes.HMM -f $VIBRANT_DB_DIR/kegg_temp.HMM $VIBRANT_DB_DIR/profile_names/VIBRANT_kegg_profiles.txt
mv $VIBRANT_DB_DIR/Pfam-A.hmm $VIBRANT_DB_DIR/Pfam-A_v32.HMM
rm -rf $VIBRANT_DB_DIR/vog_temp.HMM $VIBRANT_DB_DIR/kegg_temp.HMM $VIBRANT_DB_DIR/vog.hmm.tar.gz $VIBRANT_DB_DIR/profiles.tar.gz
hmmpress $VIBRANT_DB_DIR/VOGDB94_phage.HMM
hmmpress $VIBRANT_DB_DIR/KEGG_profiles_prokaryotes.HMM
hmmpress $VIBRANT_DB_DIR/Pfam-A_v32.HMM
#echo 'hmmfetch and hmmpressed'
chmod +x $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/scripts/*
cp -rf $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/scripts $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
chmod +x $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/VIBRANT_run.py
cp $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/VIBRANT_run.py $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
cp -r $VIBRANT_DB_DIR files $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
conda deactivate
echo '-> your VIBRANT database is now ready'
else echo '-> your VIBRANT database is ready'
fi

# WISH
WISH_DB_DIR=$database_location/wish
if [ ! -f  $database_location/.wish_db_finished ];
then
mkdir -p $WISH_DB_DIR
wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz" -P $WISH_DB_DIR
wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz" -P $WISH_DB_DIR
wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz" -P $WISH_DB_DIR
wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.1.genomic.fna.gz" -P $WISH_DB_DIR

gunzip $WISH_DB_DIR/*
cat $WISH_DB_DIR/* > $WISH_DB_DIR/viral_refseq.fna
python3 $MUDOGER_DEPENDENCIES_PATH/split-all-seq.py $WISH_DB_DIR/viral_refseq.fna $WISH_DB_DIR/viruses

for d in $WISH_DB_DIR/viruses-*fa;
do
    if grep -q phage "$d";
    then :;
else
    rm -f "$d";
fi;
done

mv $WISH_DB_DIR/viruses* $WISH_DB_DIR/phages
rm -rf $WISH_DB_DIR/viral.1.1.genomic.fna  $WISH_DB_DIR/viral.2.1.genomic.fna  $WISH_DB_DIR/viral.3.1.genomic.fna  $WISH_DB_DIR/viral.4.1.genomic.fna $WISH_DB_DIR/viral_refseq.fna
touch $database_location/.wish_db_finished
else
echo "-> your Wish database is ready"
fi


### CheckV
mkdir -p  "$database_location"/checkv
cd "$database_location"/checkv
if [ ! -d checkv-db-v1.0 ]; then
wget https://portal.nersc.gov/CheckV/checkv-db-v1.0.tar.gz
tar -zxvf checkv-db-v1.0.tar.gz
rm -fr checkv-db-v1.0.tar.gz

else echo "-> your CheckV database is ready"
fi


############################################### EUKARYOTES ###############################################
### EukCC
mkdir -p  "$database_location"/eukccdb
cd "$database_location"/eukccdb
if [ ! -d eukcc2_db_ver_1.1 ]; then
wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz
tar -xzvf eukcc2_db_ver_1.1.tar.gz
rm -fr eukcc2_db_ver_1.1.tar.gz
mv eukcc2_db_ver_1.1/* ./
rm -fr eukcc2_db_ver_1.1/

else echo "-> your EUKCC database is ready"
fi

### BUSCO eukaryotes database

mkdir -p  "$database_location"/buscodbs
cd "$database_location"/buscodbs
if [ ! -d "$database_location"/buscodbs/lineages/eukaryota_odb10 ]; then

conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/busco_env
busco --download_path "$database_location"/buscodbs --download "eukaryota"

else echo "-> your BUSCO eukaryotes database is ready"
fi
