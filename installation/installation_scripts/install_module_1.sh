################# INSTALLING MODULE 1. PRE-PROCESSING #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - METAWRAP
## - KHMER

echo '------------------------------> STARTING INSTALL MODULE 1 '
source installation/config.sh              # modified by rodolfo
source installation/installation_utils.sh  # modified by rodolfo

#echo 'stoppen sie bitte'
#echo "$PWD"
#exit 0

echo "----> mudoger dependencies path",$MUDOGER_DEPENDENCIES_ENVS_PATH


echo "### INSTALLING MODULE 1. PRE-PROCESSING ###"

## Checking if some tool already have a conda environment created
verify_if_conda_env_exist metawrap_env

### CLONE AND INSTALL METAWRAP DEPENDENCIES ###


#exit 0

echo '----> cloning metawrap github'
mkdir -p $MUDOGER_CLONED_TOOLS_PATH/metaWRAP
echo git clone $METAWRAP_GIT_URL $MUDOGER_CLONED_TOOLS_PATH/metaWRAP
git clone $METAWRAP_GIT_URL $MUDOGER_CLONED_TOOLS_PATH/metaWRAP
echo '----> done'
echo '----> creating metawrap env'
echo  create --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y python=2.7
conda create --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y python=2.7
echo '----- stoppppp'
#exit 0
echo '----> done'
echo '----> copying metawrap bin FOLDER'
echo cp -r $MUDOGER_CLONED_TOOLS_PATH/metaWRAP/bin/* $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env/bin
cp -r $MUDOGER_CLONED_TOOLS_PATH/metaWRAP/bin/* $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env/bin
echo '--> stoping now install module 1 line 38'
#exit 0
echo '----> mamba install metawrap'
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y --only-deps -c ursky -c bioconda -c conda-forge -c defaults metawrap-mg=1.3.2
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y -c bioconda fastqc
mamba install --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env -y -c bioconda trim-galore 
echo '----> done'

### INSTALLING KHMER ###
verify_if_conda_env_exist khmer_env
echo '----> creating khmer conda'
echo conda create --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/khmer_env -y python=3.6 conda-forge::r-base conda-forge::readline=6.2
conda create --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/khmer_env -y python=3.6 conda-forge::r-base conda-forge::readline=6.2
echo '----> done'
echo '----> install khmer'
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/khmer_env && pip install khmer==2.1.1 && conda deactivate    
echo '----> done'
echo '---------> END OF INSTALL MODULE 1 '

exit 0
