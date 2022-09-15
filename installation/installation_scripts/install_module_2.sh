################# MODULE 2. RECOVERY OF PROKARYOTIC MAGs #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - METAWRAP
## - GTDB-TK
## - PROKKA
## - BBTOOLS

echo "### INSTALLING MODULE 2. RECOVERY OF PROKARYOTIC MAGs ###"
source installation/config.sh            
source installation/installation_utils.sh

## Checking if some tool already have a conda environment created

################################################################
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
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y --only-deps -c ursky metawrap-mg
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y -c bioconda bwa quast megahit trim-galore fastqc metabat2
conda install -c bioconda samtools
mamba update --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y spades
echo '----> done'
fi

################################################################
## CREATE ENVIRONMENT AND INSTALLING GTDB-TK ##
verify_if_conda_env_exist gtdbtk_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/gtdbtk_env 
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/gtdbtk_env 
mamba install -y  --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/gtdbtk_env -c bioconda gtdbtk
conda deactivate
fi

################################################################
## CREATE ENVIRONMENT AND INSTALLING PROKKA ##
verify_if_conda_env_exist prokka_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/prokka_env 
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/prokka_env 
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/prokka_env -c conda-forge -c bioconda -c defaults prokka
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/prokka_env -c anaconda gawk
conda deactivate
fi
################################################################

## CREATE ENVIRONMENT AND INSTALLING GTDB-TK ##
verify_if_conda_env_exist bbtools_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/bbtools_env 
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/bbtools_env 
mamba install -y  --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/bbtools_env -c agbiome bbtools
conda deactivate
fi
################################################################
