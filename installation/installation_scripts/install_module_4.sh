################# MODULE 4. EUKARYOTIC MAGs #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - EUKREP
## - GENEMARKER-ES
## - BUSCO
## - EUKCC
source installation/config.sh              # modified by rodolfo
source installation/installation_utils.sh  # modified by rodolfo

## CREATING EUKREP ENVIRONMENT
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/eukrep_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/eukrep_env pip
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/eukrep_env -c bioconda scikit-learn==0.19.2 eukrep

## CREATING GENEMARKER-ES
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env -c bioconda perl-local-lib perl-math-utils perl-yaml perl-hash-merge perl-parallel-forkmanager perl-mce perl-threaded
#FIX HERE !!! ADD REST OF GENEMARK INSTALLATION TUTORIAL
#conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env
mkdir $MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env/GENEMARK_MAIN


## CREATING MAKER ENV
conda create y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env
#conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env anaconda::perl conda-forge::gcc conda-forge::h5py bioconda::perl-dbi bioconda::perl-module-build bioconda::perl-inline-c bioconda::perl-bit-vector bioconda::perl-io-all bioconda::perl-forks bioconda::perl-dbd-sqlite bioconda::perl-dbd-pg bioconda::perl-perl-unsafe-signals bioconda::perl-want bioconda::perl-bioperl bioconda::perl-bioperl-core bioconda::hmmer bioconda::snap bioconda::augustus bioconda::trf bioconda::blast bioconda::exonerate

#ADD REST OF MAKER 2 TUTORIAL INSTALLATION !!

## CREATING EUKCC
conda create  -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/eukcc_env -c bioconda -c biocore hmmer>3.2 pplacer python=3.7
#conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/eukcc_env -c conda-forge xorg-libxau
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/eukcc_env -c conda-forge eukcc
