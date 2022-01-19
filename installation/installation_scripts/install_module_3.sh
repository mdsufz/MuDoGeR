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
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env 
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env -c conda-forge -c bioconda virsorter#=2
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env && virsorter setup -d $MUDOGER_DEPENDENCIES_ENVS_PATH/virsorter2_env/db -j 1 && conda deactivate

## CREATING ENVIRONMENT AND INSTALLING VIRFINDER ##
verify_if_conda_env_exist virfinder_env
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virfinder_env
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/virfinder_env-c bioconda r-virfinder

## CREATING ENVIRONMENT AND INSTALLING VIBRANT ##
verify_if_conda_env_exist vibrant_env
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vcontact2_env install -y -c bioconda vibrant
conda install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env python=3 bioconda::vibrant==1.2.0 bioconda::prodigal bioconda::hmmer ostrokach::gzip conda-forge::tar conda-forge::biopython conda-forge::matplotlib anaconda::wget anaconda::pandas anaconda::seaborn anaconda::numpy anaconda::scikit-learn==0.21.3

git clone $VIBRANT_GIT_URL $MUDOGER_CLONED_TOOLS_PATH

#TODO: REVIEW AND IPROVE
conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
echo conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env


pip install pickle-mixin

VIBRANT_DB_DIR=$MUDOGER_CLONED_TOOLS_PATH/VIBRANT/databases

cho 'let us download '$VIBRANT_DB_DIR
read var

wget http://fileshare.csb.univie.ac.at/vog/vog94/vog.hmm.tar.gz -P $VIBRANT_DB_DIR
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz -P $VIBRANT_DB_DIR
wget ftp://ftp.genome.jp/pub/db/kofam/archives/2019-08-10/profiles.tar.gz -P $VIBRANT_DB_DIR

tar -xzf $VIBRANT_DB_DIR/vog.hmm.tar.gz -C $VIBRANT_DB_DIR
gunzip Pfam-A.hmm.gz -d $VIBRANT_DB_DIR
tar -xzf $VIBRANT_DB_DIR/profiles.tar.gz -C $VIBRANT_DB_DIR

for v in $VIBRANT_DB_DIR/VOG*.hmm;
do 
    cat $v >> $VIBRANT_DB_DIR/vog_temp.HMM; 
done

for k in $VIBRANT_DB_DIR/profiles/K*.hmm; 
do 
    cat $k >> $VIBRANT_DB_DIR/kegg_temp.HMM; 
done

rm -f $VIBRANT_DB_DIR/VOG0*.hmm
rm -f  $VIBRANT_DB_DIR/VOG1*.hmm
rm -f  $VIBRANT_DB_DIR/VOG2*.hmm
rm -Rf $VIBRANT_DB_DIR/profiles

hmmfetch -o $VIBRANT_DB_DIR/VOGDB94_phage.HMM -f $VIBRANT_DB_DIR/vog_temp.HMM $VIBRANT_DB_DIR/profile_names/VIBRANT_vog_profiles.txt
hmmfetch -o $VIBRANT_DB_DIR/KEGG_profiles_prokaryotes.HMM -f $VIBRANT_DB_DIR/kegg_temp.HMM $VIBRANT_DB_DIR/profile_names/VIBRANT_kegg_profiles.txt
mv $VIBRANT_DB_DIR/Pfam-A.hmm $VIBRANT_DB_DIR/Pfam-A_v32.HMM
rm -rf $VIBRANT_DB_DIR/vog_temp.HMM $VIBRANT_DB_DIR/kegg_temp.HMM $VIBRANT_DB_DIR/vog.hmm.tar.gz $VIBRANT_DB_DIR/profiles.tar.gz
hmmpress $VIBRANT_DB_DIR/VOGDB94_phage.HMM
hmmpress $VIBRANT_DB_DIR/KEGG_profiles_prokaryotes.HMM
hmmpress $VIBRANT_DB_DIR/Pfam-A_v32.HMM

chmod +x $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/scripts/*

cp -rf $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/scripts $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
chmod +x $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/VIBRANT_run.py
cp $MUDOGER_CLONED_TOOLS_PATH/VIBRANT/VIBRANT_run.py $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env
cp -r $VIBRANT_DB_DIR files $MUDOGER_DEPENDENCIES_ENVS_PATH/vibrant_env

conda deactivate

## CREATING ENVIRONMENT AND INSTALLING stampede-clustergenomes ##
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/stampede_clustergenomes_env 
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/stampede_clustergenomes_env anaconda::perl bioconda::mummer

git clone $STAMPEDE_CLUSTERGENOMES_GIT_URL $MUDOGER_CLONED_TOOLS_PATH

cp -r $MUDOGER_CLONED_TOOLS_PATH/stampede-clustergenomes/bin $MUDOGER_DEPENDENCIES_ENVS_PATH/stampede_clustergenomes_env/bin

### to run stampede-clustergenomes is necessary to activate the env and run 
### perl $your_path/mudoger_utils/cloned_tools/stampede-clustergenomes/bin/Stampede-clustergenomes.pl

## CREATING ENVIRONMENT AND INSTALLING FASTA EXTRACTION ENV
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/extract_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/extract_env anaconda::python=2.7.5

## WISH
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env conda-forge::openmp anaconda::make anaconda::cmake
git clone $WHISH_GIT_URL $MUDOGER_CLONED_TOOLS_PATH

conda activate $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env
cmake $MUDOGER_CLONED_TOOLS_PATH/WIsH
make $MUDOGER_CLONED_TOOLS_PATH/WIsH
conda deactivate

cp -r $MUDOGER_CLONED_TOOLS_PATH/WIsH $MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env/bin

WISH_DB_DIR=$MUDOGER_DEPENDENCIES_ENVS_PATH/wish_env/database

mkdir $WISH_DB_DIR


wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz" -P $WISH_DB_DIR
wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz" -P $WISH_DB_DIR
wget "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz" -P $WISH_DB_DIR

gunzip $WISH_DB_DIR/*
cat $WISH_DB_DIR/* > $WISH_DB_DIR/viral_refseq.fna
python3 $MUDOGER_DEPENDENCIES_PATH/split-all-seq.py $WISH_DB_DIR/viral_refseq.fna $WISH_DB_DIR/viruses

for d in $WISH_DB_DIR/viruses-*fa; 
do 
    if grep -q phage "$d"; 
    then :; 
else 
    rm -f "$d"; 
fi; 
done

mv $WISH_DB_DIR/viruses* $WISH_DB_DIR/phages
rm $WISH_DB_DIR/viral.1.1.genomic.fna  $WISH_DB_DIR/viral.2.1.genomic.fna  $WISH_DB_DIR/viral.3.1.genomic.fna  $WISH_DB_DIR/viral_refseq.fna

#INSTALLING CHECKV
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/checkv_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/checkv_env -c conda-forge -c bioconda checkv

#INSTALLING VCONTACT2
conda create -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vcontact2_env
mamba install -y --prefix $MUDOGER_DEPENDENCIES_ENVS_PATH/vcontact2_env python=3 pandas==0.25.1 numpy==1.16.5 -c bioconda vcontact2 mcl blast diamond prodigal

wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar -P $MUDOGER_CLONED_TOOLS_PATH

cp $MUDOGER_CLONED_TOOLS_PATH/cluster_one-1.0.jar $MUDOGER_DEPENDENCIES_ENVS_PATH/vcontact2_env/bin
