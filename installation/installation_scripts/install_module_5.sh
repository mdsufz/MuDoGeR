################# MODULE 5. BRAT ABUNDANCE #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - BOWTIE2
## - PANDASEQ
## - R
## - METAWRAP; done in initial installation
## - SAMTOOLS
## - DATAMASH
## - OTUpick
## - Prokka
## - HTSEQ

echo "### INSTALLING MODULE 5. Genome Abundance Calculation ###"

source installation/config.sh              
source installation/installation_utils.sh

############################################################################

## CREATING ENVIRONMENT, INSTALLING BRAT TOOLS ##
verify_if_conda_env_exist brat_env
if [ $PRESENT == 'yes' ]
then :;
else
mamba create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/brat_env bioconda::bowtie2 bioconda::pandaseq conda-forge::gcc conda-forge::openmpi conda-forge::r-base=3.6 bioconda::samtools bioconda::datamash

fi
############################################################################
############################################################################
## CREATING ENVIRONMENT, INSTALLING OTUpick TOOLS ##

verify_if_conda_env_exist otupick_env
if [ $PRESENT == 'yes' ]
then :;
else

conda env create --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/otupick_env --file $MUDOGER_DEPENDENCIES_PATH/otupick_dependencies/environment.yml
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/otupick_env -c anaconda gawk

chmod +x $MUDOGER_DEPENDENCIES_PATH/otupick_dependencies/*
cp -r $MUDOGER_DEPENDENCIES_PATH/otupick_dependencies/* $MUDOGER_DEPENDENCIES_ENVS_PATH/otupick_env/bin

#Checking if scripts are copied succesfully

if [ ! -s $MUDOGER_DEPENDENCIES_ENVS_PATH/otupick_env/bin/gOTUpick.sh ]; then echo "WARNING: The main script cannot be copied into conda environment. Exiting..."; exit 1; fi

if [ ! -s $MUDOGER_DEPENDENCIES_ENVS_PATH/otupick_env/bin/organize-bins-tax.py ]; then echo "WARNING: The auxiliary script organize-bins-tax.py cannot be copied into conda environment. Exiting..."; exit 1; fi

if [ ! -s $MUDOGER_DEPENDENCIES_ENVS_PATH/otupick_env/bin/aniSplitter.R ]; then echo "WARNING: The auxiliary script aniSplitter.R cannot be copied into conda environment. Exiting..."; exit 1; fi

if [ ! -s $MUDOGER_DEPENDENCIES_ENVS_PATH/otupick_env/bin/pick_rep.R ]; then echo "WARNING: The auxiliary script pick_rep.R cannot be copied into conda environment. Exiting..."; exit 1; fi

if [ ! -s $MUDOGER_DEPENDENCIES_ENVS_PATH/otupick_env/bin/summarize-anisplitter-results.py ]; then echo "WARNING: The auxiliary script summarize-anisplitter-results.py cannot be copied into conda environment. Exiting..."; exit 1; fi

fi
############################################################################

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

############################################################################

## CREATE ENVIRONMENT AND INSTALLING HTSEQ ##
verify_if_conda_env_exist htseq_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/htseq_env python=3.7
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/htseq_env
pip install HTSeq
conda deactivate
fi


############################################################################

## CREATE ENVIRONMENT AND INSTALLING Cov and TPM calculation ##
verify_if_conda_env_exist cov_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/cov_env python=3.7
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/cov_env
pip install pandas
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/cov_env -c anaconda gawk
conda deactivate
fi
