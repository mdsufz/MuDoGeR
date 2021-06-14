# this scripts asks for the desired location of installation for the databases.
# given the user input, the config file will be edited and all databases via wget, curl, etc and zipping


source config


#while: do read $input


### CheckM
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf *.tar.gz
rm *.gz
checkm data setRoot     # CheckM will prompt to to chose your storage location
# On newer versions of CheckM, you would run:
checkm data setRoot /path/to/your/dir/MY_CHECKM_FOLDER


### GTDB-tk
# download database
cd 
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz -P 
gunzip -xvf /gtdbtk_r95_data.tar.gz
