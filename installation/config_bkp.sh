"$( echo "$MUDOGER_CONDA_ENVIRONMENT_PATH")"
"$( echo MUDOGER_CONDA_ENVIRONMENT_PATH)"
"$MUDOGER_CONDA_ENVIRONMENT_PATH"
$MUDOGER_CONDA_ENVIRONMENT_PATH
#echo "------> loading config.sh"
#source $(dirname $0)/installation_utils.sh

#conda activate
#BASE_CONDA=$CONDA_PREFIX
#MUDOGER_CONDA_ENVIRONMENT_PATH=$CONDA_PREFIX/envs/mudoger_env
#conda deactivate

MUDOGER_DEPENDENCIES_PATH=$MUDOGER_CONDA_ENVIRONMENT_PATH/dependencies

MUDOGER_DEPENDENCIES_ENVS_PATH=$MUDOGER_DEPENDENCIES_PATH/conda/envs
MUDOGER_CLONED_TOOLS_PATH=$MUDOGER_DEPENDENCIES_PATH/cloned_tools
MUDOGER_INSTALLATION_SCRIPTS_PATH=$MUDOGER_DEPENDENCIES_PATH/installation_scripts

METAWRAP_GIT_URL=https://github.com/bxlab/metaWRAP.git
UBIN_GIT_URL=https://github.com/ProbstLab/uBin-helperscripts.git
VIBRANT_GIT_URL=https://github.com/AnantharamanLab/VIBRANT
STAMPEDE_CLUSTERGENOMES_GIT_URL=https://bitbucket.org/MAVERICLab/stampede-clustergenomes.git
WHISH_GIT_URL=https://github.com/soedinglab/WIsH.git

DEPENDENCIES_SCRIPTS_PATH=installation/dependencies
DEPENDENCIES_SCRIPTS_PATH=installation/dependencies # modified by rodolfo

INSTALLATION_SCRIPTS_PATH=installation/installation_scripts
