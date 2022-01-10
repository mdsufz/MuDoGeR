################# MODULE 2. RECOVERY OF PROKARYOTIC MAGs #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - GTDB-TK
## - PROKKA
## - UBIN

echo "### INSTALLING MODULE 2. RECOVERY OF PROKARYOTIC MAGs ###"
source installation/config.sh              # modified by rodolfo
source installation/installation_utils.sh  # modified by rodolfo
## Checking if some tool already have a conda environment created
verify_if_conda_env_exist gtdbtk_env
verify_if_conda_env_exist prokka_env
verify_if_conda_env_exist ubin_env

## CREATE ENVIRONMENT AND INSTALLING GTDB-TK ##
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/gtdbtk_env 
mamba install -y  --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/gtdbtk_env -c bioconda gtdbtk

## CREATE ENVIRONMENT AND INSTALLING PROKKA ##
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/prokka_env 
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/prokka_env  -c conda-forge -c bioconda -c defaults prokka

## CLONE UBIN AND CREATING ENVIRONMENT BASED ON .YAML FILE ##
#git clone $UBIN_GIT_URL $MUDOGER_CLONED_TOOLS_PATH
#conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/ubin_env
#mamba install -y  --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/ubin_env -f $MUDOGER_CLONED_TOOLS_PATH/uBin-helperscripts/uBin_wrapper_reqs.yaml