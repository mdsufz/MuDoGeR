# this scripts asks for the desired location of installation for the databases.
# given the user input, the config file will be edited and all databases via wget, curl, etc and zipping


source config


while: do read $input


### GTDB-tk
# download database
cd 
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz -P 
gunzip -xvf /gtdbtk_r95_data.tar.gz
