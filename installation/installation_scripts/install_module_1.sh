################# INSTALLING MODULE 1. PRE-PROCESSING #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - METAWRAP
## - KHMER

echo '------------------------------> STARTING INSTALL MODULE 1 '
source installation/config.sh              
source installation/installation_utils.sh  


echo "----> mudoger dependencies path",$MUDOGER_DEPENDENCIES_ENVS_PATH


echo "### INSTALLING MODULE 1. PRE-PROCESSING ###"

## Checking if some tool already have a conda environment created
### CLONE AND INSTALL METAWRAP DEPENDENCIES ###
verify_if_conda_env_exist metawrap_env
if [ $PRESENT == 'yes' ]
then :;
else
echo '----> cloning metawrap github'
mkdir -p $MUDOGER_CLONED_TOOLS_PATH/metaWRAP
echo git clone $METAWRAP_GIT_URL $MUDOGER_CLONED_TOOLS_PATH/metaWRAP
git clone $METAWRAP_GIT_URL $MUDOGER_CLONED_TOOLS_PATH/metaWRAP
echo '----> done'
echo '----> creating metawrap env'
echo  create --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y python=2.7
conda create --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y python=2.7

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels ursky
echo '----> done'

echo '----> copying metawrap bin FOLDER'
echo cp -r $MUDOGER_CLONED_TOOLS_PATH/metaWRAP/bin/* $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env/bin
cp -r $MUDOGER_CLONED_TOOLS_PATH/metaWRAP/bin/* $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env/bin

echo '----> mamba install metawrap'
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y --only-deps -c ursky metawrap-mg=1.3.2
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y -c bioconda samtools=1.9 --force-reinstall
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y -c bioconda bwa quast megahit trim-galore fastqc metabat2 spades
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y concoct=1.0
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y checkm-genome
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y maxbin2
mamba update --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y spades
echo '----> done'
fi

### INSTALLING KHMER ###
verify_if_conda_env_exist khmer_env
if [ $PRESENT == 'yes' ]
then :;
else
echo '----> creating khmer env'
echo conda create --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/khmer_env -y python=3.8 conda-forge::r-base conda-forge::readline=6.2
conda create --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/khmer_env -y python=3.8 conda-forge::r-base conda-forge::readline=6.2
echo '----> done'
echo '----> install khmer'
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/khmer_env && pip install khmer==2.1.1 && conda deactivate    
echo '----> done'
echo '---------> END OF INSTALL MODULE 1 '
fi


