################# MODULE 2. RECOVERY OF PROKARYOTIC MAGs #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - GTDB-TK
## - PROKKA
## - UBIN
## - BBTOOLS

echo "### INSTALLING MODULE 2. RECOVERY OF PROKARYOTIC MAGs ###"
source installation/config.sh              # modified by rodolfo
source installation/installation_utils.sh  # modified by rodolfo
## Checking if some tool already have a conda environment created


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

## CLONE UBIN AND CREATING ENVIRONMENT BASED ON .YAML FILE ##
#git clone $UBIN_GIT_URL $MUDOGER_CLONED_TOOLS_PATH
#conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/ubin_env
#mamba install -y  --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/ubin_env -f $MUDOGER_CLONED_TOOLS_PATH/uBin-helperscripts/uBin_wrapper_reqs.yaml
