################# MODULE 5. BRAT ABUNDANCE #################

### TOOLS THAT WILL BE INSTALLED IN THIS MODULE ###
## - BOWTIE2
## - PANDASEQ
## - R
## - METAWRAP; done in initial installation
## - SAMTOOLS
## - DATAMASH

echo "### INSTALLING MODULE 5. Genome Abundance Calculation ###"

source installation/config.sh              
source installation/installation_utils.sh

############################################################################

## CREATING ENVIRONMENT, INSTALLING BRAT TOOLS ##
verify_if_conda_env_exist brat_env
if [ $PRESENT == 'yes' ]
then :;
else
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/brat_env bioconda::bowtie2 bioconda::pandaseq conda-forge::gcc conda-forge::openmpi conda-forge::r-base=3.6 bioconda::samtools bioconda::datamash

fi
############################################################################
