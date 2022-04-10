######### COMPLETING INSTALLATION FOR GENEMARK #########
#0) Activate mudoger_env
conda activate mudoger_env
config_file="$(which config.sh)"
source "$config_file"

MUDOGER_DEPENDENCIES_ENVS_PATH=$MUDOGER_DEPENDENCIES_ENVS_PATH
# 1) ACCESS GENEMARK WEBPAGE
http://opal.biology.gatech.edu/GeneMark/license_download.cgi

# 2) SELECT OPTIONS
# GeneMark-ES/ET/EP ver 4.68_lic and LINUX 64

# 3) FILL IN THE CREDENTIALS WITH YOUR NAME, E-MAIL, INSTITUTION, ETC... 

# 4) CLICK ON 'I agree the terms of this license agreement'

# 5) DOWNLOAD THE 64_bit key and the program files 
#ex_files: wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_HZzc0/gmes_linux_64.tar.gz
#ex_key: wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_HZzc0/gm_key_64.gz

# 5.1) You should have the following files:
#gmes_linux_64.tar.gz  gm_key_64.gz

# 6) RETRIEVE PATH TO GENEMARK ENV (this environment is already created and ready with install_module_4.sh)
genemarker_environment="$MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env"

# 7) PLACE BOTH FILES INSIDE A <CONDA_GENEMARK>/GENEMARK_MAIN. 
cp /path/to/gmes_linux_64.tar.gz  $genemarker_environment/GENEMARK_MAIN
cp /path/to/gm_key_64.gz  $genemarker_environment/GENEMARK_MAIN
# IT SHOULD LOOK LIKE THIS:
ls  $genemarker_environment/GENEMARK_MAIN
gmes_linux_64.tar.gz  gm_key_64.gz


# 8) ENTER INSIDE FOLDER AND DECOMPRESS FILES
cd  $genemarker_environment/GENEMARK_MAIN
gunzip *

# 9) COPY AND RENAME KEY FILE TO HOME FOLDER 
cp gm_key_64 ~/.gm_key

# 10) DECOMPRESS TOOL 
tar -xvf gmes_linux_64.tar 

# 11) ACTIVATE THE genemark_env CONDA 
conda activate $genemarker_environment

# 12) GET PERL PATH
perl_path="$(which perl)"

cd - #return
# 13) CONFIGURE GENEMARK TO USE THE CORRECT PERL
cd $genemarker_environment/GENEMARK_MAIN/gmes_linux_64
perl change_path_in_perl_scripts.pl "$perl_path"

cd - #return
# 14) Return to mudoger_env
conda deactivate

######### COMPLETING INSTALLATION FOR MAKER2 #########

# 1) DOWNLOAD MAKER (Registration required)

#Go to http://weatherby.genetics.utah.edu/cgi-bin/registration/maker_license.cgi and fill out your information to download the tool
#Copy the download link provided (Usually by Right-click -> 'Copy Link') and type: 

wget *Download link*
#ex: wget http://weatherby.genetics.utah.edu/maker_downloads/1BF1/66D1/3FAA/FF20FB345D6221AA8A338D1B9D8A/maker-3.01.04.tgz

# 2) RETRIEVE PATH TO MUDOGER_ENV (this environment is already created and ready with install.sh)
maker2_environment="$MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env" 

# 3) PLACE THE DOWNLOADED MAKER FILE INSIDE <CONDA_MAKER2>/MAKER_MAIN. 
cp maker-3.01.04.tgz $maker2_environment/MAKER_MAIN

# 4) ENTER INSIDE FOLDER AND DECOMPRESS FILES
cd  $maker2_environment/MAKER_MAIN
tar xvfz maker-3.01.04.tgz

# 5) CREATE A EXE FOLDER INSIDE THE MAKER EXTRACTED FOLDER
mkdir $maker2_environment/MAKER_MAIN/maker/exe

# 6) DOWNLOAD DFAM DATABASE. You should have enough space for this database where you installed the conda envs (~97GB)

wget -P $maker2_environment/MAKER_MAIN/ https://www.dfam.org/releases/Dfam_3.5/families/Dfam.h5.gz
gunzip Dfam.h5.gz

# 7) DOWNLOAD AND MOVE SNAP

wget -P $maker2_environment/MAKER_MAIN/ http://korflab.ucdavis.edu/Software/snap-2013-11-29.tar.gz
tar xvfz snap-2013-11-29.tar.gz
mv snap $maker2_environment/MAKER_MAIN/maker/exe

# 8) COPY SOME SCRIPTS FROM YOUR GENEMARK INSTALLATION (Installation done previously, just copied the necessaries bin)

cp -r $MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env/GENEMARK_MAIN/gmes_linux_64/* $MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env/bin/

# 9) DOWNLOAD REPEATMASKER

wget -P $maker2_environment/MAKER_MAIN/ https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz

# 10) EXTRACT REPEATMASKER

tar xvfz RepeatMasker-4.1.2-p1.tar.gz

# 11) MOVE DFAM DATABASE TO LIBRARIES FOLDER IN REPEATMASKER

mv -f $maker2_environment/MAKER_MAIN/Dfam.h5 $maker2_environment/MAKER_MAIN/RepeatMasker/Libraries/ #(Remember to force yes)

# 12) MOVE REPEATMASKER TO EXE FOLDER IN MAKER

mv $maker2_environment/MAKER_MAIN/RepeatMasker $maker2_environment/MAKER_MAIN/maker/exe/

cd - #return
# 13) INSTALL REPEATMASKER

conda activate $maker2_environment

cd $maker2_environment/MAKER_MAIN/maker/exe/RepeatMasker/
perl ./configure 

#Manually give the paths 
#Ideally it should automatically find everything. If that is the case, just hit enter when prompt
#Hit enter if TRF path is correct
#Select 3 (HMMER) as search engine and set as default (Y)
#Hit 5 to be done and wait databse configuration

cd - #return
# 14) CONFIGURE MAKER

cd $maker2_environment/MAKER_MAIN/maker/src

perl Build.PL #(Currently answer N to set MPI)
##If any Perl dependency was not properly installed, for any reason, please try the following

cpanm FAILED::DEPENDENCY #ex: cpanm Bit::Vector

# 15) INSTALL MAKER
./Build install

export ZOE="$maker2_environment/MAKER_MAIN/maker/exe/snap/Zoe"
export AUGUSTUS_CONFIG_PATH="$maker2_environment/config"

# 16) ADD MAKER TO PATH

export PATH=$maker2_environment/MAKER_MAIN/maker/bin:$PATH

cd - #return
conda deactivate
