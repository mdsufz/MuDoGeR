################# MODULE 4. EUKARYOTIC MAGs #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - EUKREP
## - GENEMARKER-ES
## - MAKER2
## - EUKCC
## - BUSCO

echo "### INSTALLING MODULE 4. RECOVERY OF EUKARYOTIC BINS ###"

source installation/config.sh              
source installation/installation_utils.sh

############################################################################

## CREATING ENVIRONMENT, INSTALLING EUKREP ##
verify_if_conda_env_exist eukrep_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/eukrep_env -c bioconda scikit-learn==0.19.2 eukrep

fi
############################################################################

############################################################################

## CREATING ENVIRONMENT, INSTALLING GENEMARKER-ES ##
verify_if_conda_env_exist genemarker_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env
#mamba have problems with perl
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env -c bioconda perl-local-lib perl-math-utils perl-yaml perl-hash-merge perl-parallel-forkmanager perl-mce perl-threaded
mkdir $MUDOGER_DEPENDENCIES_ENVS_PATH/genemarker_env/GENEMARK_MAIN

fi
############################################################################

############################################################################

## CREATING ENVIRONMENT, INSTALLING MAKER2 ##
verify_if_conda_env_exist maker2_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env
#mamba has problems with perl
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env anaconda::perl conda-forge::gcc conda-forge::h5py bioconda::perl-dbi bioconda::perl-module-build bioconda::perl-inline-c bioconda::perl-bit-vector bioconda::perl-io-all bioconda::perl-forks bioconda::perl-dbd-sqlite bioconda::perl-dbd-pg bioconda::perl-perl-unsafe-signals bioconda::perl-want bioconda::perl-bioperl bioconda::perl-bioperl-core bioconda::hmmer bioconda::snap bioconda::augustus bioconda::trf bioconda::blast bioconda::exonerate
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env
cpanm List::MoreUtils
mkdir $MUDOGER_DEPENDENCIES_ENVS_PATH/maker2_env/MAKER_MAIN
conda deactivate

fi
############################################################################

############################################################################

## CREATING ENVIRONMENT, INSTALLING EUKCC ##
verify_if_conda_env_exist eukcc_env
if [ $PRESENT == 'yes' ]
then :;
else

conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/eukcc_env -c conda-forge -c bioconda "eukcc>=2"

fi

############################################################################

## CREATING ENVIRONMENT, INSTALLING BUSCO ##
verify_if_conda_env_exist busco_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/busco_env -c conda-forge -c bioconda busco=5.3.0

fi
############################################################################
