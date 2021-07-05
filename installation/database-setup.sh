# this scripts asks for the desired location of installation for the databases.
# given the user input, the config file will be edited and all databases via wget, curl, etc and zipping

database_location="$1"

mkdir "$database_location"

source config



### CheckM
mkdir  "$database_location"/"checkm"
cd  "$database_location"/"checkm"
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf *.tar.gz
rm *.gz
checkm data setRoot     # CheckM will prompt to to chose your storage location
# On newer versions of CheckM, you would run:
checkm data setRoot /path/to/your/dir/MY_CHECKM_FOLDER


### GTDB-tk
# download database
cd $database_location
mkdir  "$database_location"/"gtdbtk"
cd "$database_location"/"gtdbtk"
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz -P 
gunzip -xvf /gtdbtk_r95_data.tar.gz
echo  GTDBTK_DATA_PATH=$database_location/gtdbtk/gtdbtk_r95_data >> config
