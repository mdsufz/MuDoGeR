# this scripts asks for the desired location of installation for the databases.
# given the user input, the config file will be edited and all databases via wget, curl, etc and zipping

database_location="$1"

mkdir "$database_location"

source config
############################################### PROKARYOTES ###############################################
### CheckM
mkdir -p "$database_location"/checkm
cd  "$database_location"/checkm
if [ ! -d checkm_data_2015_01_16 ]; then
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf checkm_data_2015_01_16.tar.gz
rm checkm_data_2015_01_16.tar.gz
# On newer versions of CheckM, you would run:
#checkm data setRoot /path/to/your/dir/$MY_CHECKM_FOLDER
CHECKM_DB="$database_location"/"checkm"/checkm_data_2015_01_16
echo  CHECKM_DB="$CHECKM_DB" >> config
else echo "-> your CheckM database is ready"


### GTDB-tk
mkdir -p  "$database_location"/"gtdbtk"
cd "$database_location"/"gtdbtk"
if [ ! -d gtdbtk_r95_data ]; then
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz -P 
gunzip -xvf gtdbtk_r95_data.tar.gz
rm gtdbtk_r95_data.tar.gz
echo  GTDBTK_DATA_PATH=$database_location/gtdbtk/gtdbtk_r95_data >> config
else echo "-> your GTDBtk database is ready"

############################################### VIRUSES ###############################################
### CheckV
mkdir -p  "$database_location"/checkv
cd "$database_location"/"gtdbtk"
wget https://portal.nersc.gov/CheckV/checkv-db-v1.0.tar.gz
tar -zxvf checkv-db-v1.0.tar.gz
