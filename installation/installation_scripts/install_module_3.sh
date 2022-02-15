################# MODULE 3. UNCULTIVATED VIRAL MAGs #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - VIRSORTER 2
## - VIRFINDER
## - VIBRANT
## - stampede-clustergenomes
## - WIsH
## - CHECKV
## - VCONTACT2
## - ClusterONE

echo "### INSTALLING MODULE 3. RECOVERY OF UVIGS ###"
source installation/config.sh              # modified by rodolfo
source installation/installation_utils.sh  # modified by rodolfo
## Checking if some tool already have a conda environment created


############################################################################
## CREATING ENVIRONMENT, INSTALLING VIRSORTER 2 AND SETUP DATABASE ##
verify_if_conda_env_exist virsorter2_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env -c conda-forge -c bioconda virsorter#=2
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env && virsorter setup -d $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env/db -j 1 && conda deactivate
fi
############################################################################
## CREATING ENVIRONMENT AND INSTALLING VIRFINDER ##
verify_if_conda_env_exist virfinder_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virfinder_env
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virfinder_env -c bioconda r-virfinder
fi

############################################################################
## CREATING ENVIRONMENT AND INSTALLING VIBRANT ##
verify_if_conda_env_exist vibrant_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vcontact2_env install -y -c bioconda vibrant
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env python=3 bioconda::vibrant==1.2.0 bioconda::prodigal bioconda::hmmer ostrokach::gzip conda-forge::tar conda-forge::biopython conda-forge::matplotlib anaconda::wget anaconda::pandas anaconda::seaborn anaconda::numpy anaconda::scikit-learn==0.21.3
git clone $VIBRANT_GIT_URL $MUDOGER_CLONED_TOOLS_PATH
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
echo conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
pip install pickle-mixin
fi

############################################################################
## CREATING ENVIRONMENT AND INSTALLING stampede-clustergenomes ##
verify_if_conda_env_exist stampede_clustergenomes_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/stampede_clustergenomes_env 
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/stampede_clustergenomes_env anaconda::perl bioconda::mummer
git clone $STAMPEDE_CLUSTERGENOMES_GIT_URL $MUDOGER_CLONED_TOOLS_PATH/stampede-clustergenomes
cp -rf $MUDOGER_CLONED_TOOLS_PATH/stampede-clustergenomes/bin/* $MUDOGER_DEPENDENCIES_ENVS_PATH/stampede_clustergenomes_env/bin
fi

############################################################################
## CREATING ENVIRONMENT AND INSTALLING FASTA EXTRACTION ENV
verify_if_conda_env_exist extract_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/extract_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/extract_env python=2
fi


#verify_if_conda_env_exist extract_env
#if [ $PRESENT == 'yes' ]
#then :;
#else
#conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env
#mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env conda-forge::openmp anaconda::make anaconda::cmake
#fi

############################################################################
## WISH
verify_if_conda_env_exist wish_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env conda-forge::openmp anaconda::make anaconda::cmake
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env
git clone $WISH_GIT_URL $MUDOGER_CLONED_TOOLS_PATH/WIsH
cd $MUDOGER_CLONED_TOOLS_PATH/WIsH
cmake .
make
cp WIsH $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env/bin
cd -
conda deactivate
WISH_DB_DIR=$MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env/database
mkdir $WISH_DB_DIR
fi


############################################################################
#INSTALLING CHECKV
verify_if_conda_env_exist checkv_env
if [ $PRESENT == 'yes' ]
then ;
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/checkv_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/checkv_env -c conda-forge -c bioconda checkv
else :
fi

############################################################################
#INSTALLING VCONTACT2
verify_if_conda_env_exist vcontact2_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vcontact2_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vcontact2_env  -c conda-forge python=3 pandas==0.25.1 numpy==1.16.5 
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vcontact2_env  -c bioconda vcontact2 mcl blast diamond prodigal
wget --no-check http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -P $MUDOGER_CLONED_TOOLS_PATH
cp $MUDOGER_CLONED_TOOLS_PATH/cluster_one-1.0.jar $MUDOGER_DEPENDENCIES_ENVS_PATH/vcontact2_env/bin
fi
