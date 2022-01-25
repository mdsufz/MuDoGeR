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



## CREATING ENVIRONMENT, INSTALLING VIRSORTER 2 AND SETUP DATABASE ##
verify_if_conda_env_exist virsorter2_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env -c conda-forge -c bioconda virsorter#=2
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env && virsorter setup -d $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env/db -j 1 && conda deactivate
fi

## CREATING ENVIRONMENT AND INSTALLING VIRFINDER ##
verify_if_conda_env_exist virfinder_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virfinder_env
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virfinder_env-c bioconda r-virfinder
fi

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


### to run stampede-clustergenomes is necessary to activate the env and run 
### perl $your_path/mudoger_utils/cloned_tools/stampede-clustergenomes/bin/Stampede-clustergenomes.pl

## CREATING ENVIRONMENT AND INSTALLING FASTA EXTRACTION ENV
verify_if_conda_env_exist extract_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/extract_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/extract_env anaconda::python=2.7.5
fi

## WISH
verify_if_conda_env_exist extract_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env conda-forge::openmp anaconda::make anaconda::cmake

$WHISH_GIT_URL $MUDOGER_CLONED_TOOLS_PATH
fi

verify_if_conda_env_exist wish_env
if [ $PRESENT == 'yes' ]
then :;
else
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env
cmake $MUDOGER_CLONED_TOOLS_PATH/WIsH
make $MUDOGER_CLONED_TOOLS_PATH/WIsH
conda deactivate
cp -r $MUDOGER_CLONED_TOOLS_PATH/WIsH $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env/bin
WISH_DB_DIR=$MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env/database
mkdir $WISH_DB_DIR
fi

#if [
#wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz" -P $WISH_DB_DIR
#wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz" -P $WISH_DB_DIR
#wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz" -P $WISH_DB_DIR
#gunzip $WISH_DB_DIR/*
#cat $WISH_DB_DIR/* > $WISH_DB_DIR/viral_refseq.fna
#python3 $MUDOGER_DEPENDENCIES_PATH/split-all-seq.py $WISH_DB_DIR/viral_refseq.fna $WISH_DB_DIR/viruses

#for d in $WISH_DB_DIR/viruses-*fa; 
#do 
#    if grep -q phage "$d"; 
#    then :; 
#else 
#    rm -f "$d"; 
#fi; 
#done

#mv $WISH_DB_DIR/viruses* $WISH_DB_DIR/phages
#rm -rf $WISH_DB_DIR/viral.1.1.genomic.fna  $WISH_DB_DIR/viral.2.1.genomic.fna  $WISH_DB_DIR/viral.3.1.genomic.fna  $WISH_DB_DIR/viral_refseq.fna

#INSTALLING CHECKV
verify_if_conda_env_exist checkv_env
if [ $PRESENT == 'yes' ]
then :;
else :
#conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/checkv_env
#mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/checkv_env -c conda-forge -c bioconda checkv
fi

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
