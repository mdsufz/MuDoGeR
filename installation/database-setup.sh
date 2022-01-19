#!/bin/bash

# this scripts asks for the desired location of installation for the databases.
# given the user input, the config file will be edited and all databases via wget, curl, etc and zipping

database_location="$1"
#config_file=MuDoGeR/installation/config #Path problem! The script does not know where this config is


#touch "$config_file"

mkdir "$database_location"
conda activate mudoger_env
config_file="$(which config.sh)"
source "$config_file"
############################################### PROKARYOTES ###############################################
### CheckM
echo 'installing checkm database ...'
mkdir -p "$database_location"/checkm
cd  "$database_location"/checkm
if [ ! -f selected_marker_sets.tsv ]; then
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf checkm_data_2015_01_16.tar.gz
rm -fr checkm_data_2015_01_16.tar.gz
# On newer versions of CheckM, you would run:
#checkm data setRoot /path/to/your/dir/$MY_CHECKM_FOLDER
CHECKM_DB="$database_location"/checkm #Fixed? we need to test
echo CHECKM_DB="$CHECKM_DB" >> "$config_path"
else echo "-> your CheckM database is ready"
fi


### GTDB-tk
mkdir -p  "$database_location"/"gtdbtk"
cd "$database_location"/"gtdbtk"
if [ ! -d release* ]; then
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
rm -fr gtdbtk_data.tar.gz
#echo  GTDBTK_DATA_PATH="$database_location"/gtdbtk/gtdbtk_r95_data >> "$config_file" ##FIX HERE
echo GTDBTK_DATA_PATH="$(ls "$database_location"/gtdbtk/release*)" >>  "$config_file" # fixed? we need to test
else echo "-> your GTDBtk database is ready"
fi

############################################### VIRUSES ###############################################
### VIBRANT
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
echo conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env

pip install pickle-mixin

VIBRANT_DB_DIR=$MUDOGER_CLONED_TOOLS_PATH/VIBRANT/databases

cho 'let us download '$VIBRANT_DB_DIR
read var

wget http://fileshare.csb.univie.ac.at/vog/vog94/vog.hmm.tar.gz -P $VIBRANT_DB_DIR
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz -P $VIBRANT_DB_DIR
wget ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/profiles.tar.gz -P $VIBRANT_DB_DIR

tar -xzf $VIBRANT_DB_DIR/vog.hmm.tar.gz -C $VIBRANT_DB_DIR
gunzip Pfam-A.hmm.gz -d $VIBRANT_DB_DIR
tar -xzf $VIBRANT_DB_DIR/profiles.tar.gz -C $VIBRANT_DB_DIR

for v in $VIBRANT_DB_DIR/VOG*.hmm;
do 
    cat $v >> $VIBRANT_DB_DIR/vog_temp.HMM; 
done

for k in $VIBRANT_DB_DIR/profiles/K*.hmm; 
do 
    cat $k >> $VIBRANT_DB_DIR/kegg_temp.HMM; 
done

rm -f $VIBRANT_DB_DIR/VOG0*.hmm
rm -f  $VIBRANT_DB_DIR/VOG1*.hmm
rm -f  $VIBRANT_DB_DIR/VOG2*.hmm
rm -Rf $VIBRANT_DB_DIR/profiles

hmmfetch -o $VIBRANT_DB_DIR/VOGDB94_phage.HMM -f $VIBRANT_DB_DIR/vog_temp.HMM $VIBRANT_DB_DIR/profile_names/VIBRANT_vog_profiles.txt
hmmfetch -o $VIBRANT_DB_DIR/KEGG_profiles_prokaryotes.HMM -f $VIBRANT_DB_DIR/kegg_temp.HMM $VIBRANT_DB_DIR/profile_names/VIBRANT_kegg_profiles.txt
mv $VIBRANT_DB_DIR/Pfam-A.hmm $VIBRANT_DB_DIR/Pfam-A_v32.HMM
rm -rf $VIBRANT_DB_DIR/vog_temp.HMM $VIBRANT_DB_DIR/kegg_temp.HMM $VIBRANT_DB_DIR/vog.hmm.tar.gz $VIBRANT_DB_DIR/profiles.tar.gz
hmmpress $VIBRANT_DB_DIR/VOGDB94_phage.HMM
hmmpress $VIBRANT_DB_DIR/KEGG_profiles_prokaryotes.HMM
hmmpress $VIBRANT_DB_DIR/Pfam-A_v32.HMM

chmod +x $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/scripts/*

cp -rf $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/scripts $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
chmod +x $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/VIBRANT_run.py
cp $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/VIBRANT_run.py $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
cp -r $VIBRANT_DB_DIR files $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env

conda deactivate




### CheckV
mkdir -p  "$database_location"/checkv
cd "$database_location"/checkv
if [ ! -d checkv-db-v1.0 ]; then
wget https://portal.nersc.gov/CheckV/checkv-db-v1.0.tar.gz
tar -zxvf checkv-db-v1.0.tar.gz
rm -fr checkv-db-v1.0.tar.gz
#ADD CHECKV DATABASE PATH TO CONFIG FILE
CHECKVDB="$database_location"/checkv/checkv-db-v1.0
echo CHECKVDB="$CHECKVDB" >> "$config_file"
else echo "-> your CheckV database is ready"
fi

############################################### EUKARYOTES ###############################################




