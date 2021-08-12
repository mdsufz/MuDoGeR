#!/bin/bash

## Try to go the cloned_tools directory
## If cloned_tools folder does not exist, will
## create it and enter into this new directory
go_to_cloned_tools_folder() {
	cd cloned_tools || mkdir cloned_tools && cd cloned_tools
}

## Verify if a given conda environment exist
## Assign the return to a variable called result
## If result=0: EXIST
## If result=1: DOES NOT EXIST
verify_if_conda_env_exist () {
	local conda_env_name=$1
	local conda_list=$(conda list -n $conda_env_name)
	if [ "$conda_list" ];
	then
		result=0
	else
		result=1
	fi
}

## Create a Metawrap conda environment
clone_and_install_metawrap() {
	#INSTALLING METAWRAP
	go_to_cloned_tools_folder
	git clone https://github.com/bxlab/metaWRAP.git
	#->configure the config-metawrap file to pointo to the database
	conda create -y -n metawrap-env python=2.7
	cd metaWRAP/bin
	cp -r . $conda_path/envs/metawrap-env/bin
	cd ../../..
	conda activate metawrap-env
	#conda config --add channels defaults
	#conda config --add channels conda-forge
	#conda config --add channels bioconda
	#conda config --add channels ursky
	conda install -y --only-deps -c ursky metawrap-mg=1.3.2
	metawrap
	conda deactivate
}

#making a directory for mudoger's utilities and change directory to it
mkdir mudoger_utils
cd mudoger_utils

#creating conda path variable
echo "Please, enter with your miniconda path (ex. \"/home/profile_name/miniconda3\")"
read conda_path
echo "your conda path: $conda_path"

echo "Do you want to install all the required tools for the Module 1 (Pre-processing)?[y/n]"
while :
do
	read choose
	if [ $choose = y -o $choose = Y ];
	then
		#INSTALLING METAWRAP
		clone_and_install_metawrap
		
		## INSTALLING KHMER
		conda create -y -n khmer-env python=3.6
		conda activate khmer-env
		pip install khmer==2.1.1
		conda install -y -c conda-forge r-base
		conda deactivate
		break
	elif [ $choose = n -o $choose = N ]
	then
		echo "installation denied..."
		break
	else
		echo "command not found, please try again"
	fi
done

echo "Do you want to install all the required tools for the Module 2 (Recovery of Procaryotic Metagenome-Assembled Genomes)?[y/n]"
while :
do
	read choose
	if [ $choose = y -o $choose = Y ];
	then	

		#INSTALLING METAWRAP IF NECESSARY
		verify_if_conda_env_exist "metawrap-env"
		if [ "$result" = 1 ];
		then
			clone_and_install_metawrap
		fi

		#INSTALLING GTDB-TK
		conda create -y -n gtdbtk-env
		conda install -n gtdbtk-env -y -c bioconda gtdbtk

		#INSTALLING UBIN
		go_to_cloned_tools_folder
		git clone https://github.com/ProbstLab/uBin-helperscripts.git
		cd uBin-helperscripts
		conda env create -n ubin-env -f uBin_wrapper_reqs.yaml
		cd ../..

		break
	elif [ $choose = n -o $choose = N ]
	then
		echo "installation denied..."
		break
	else
		echo "command not found, please try again"
	fi
done

echo "Do you want to install all the required tools for the Module 3 (Recovery of Uncutivated Viral Genomes)?[y/n]"
while :
do
	read choose
	if [ $choose = y -o $choose = Y ];
	then
		#INSTALLING VIRSORTER 1
		#conda create --name virsorter-env -y -c bioconda mcl=14.137 muscle blast perl-bioperl perl-file-which hmmer=3.1b2 perl-parallel-forkmanager perl-list-moreutils diamond=0.9.14
		#go_to_cloned_tools_folder
		#git clone https://github.com/simroux/VirSorter.git
		#cd VirSorter/Scripts
		#conda install -n virsorter-env -y -c anaconda make
		#conda activate virsorter-env
		#make clean
		#make
		#conda deactivate
		#cp -r . $conda_path/envs/virsorter-env/bin
		#cd ..
		#cp wrapper_phage_contigs_sorter_iPlant.pl $conda_path/envs/virsorter-env/bin
		#conda install -n virsorter-env -y -c bioconda metagene_annotator
		#cd ../..
		
		#INSTALLING VIRSORTER 2
		conda create -n virsorter2-env -y -c conda-forge -c bioconda virsorter=2
		conda activate virsorter2-env
		conda deactivate


		#INSTALLING VIRFINDER
		#conda create -y -n virfinder-env
		#conda install -n virfinder-env -y -c bioconda r-virfinder

		#INSTALLING VIBRANT
		#conda create -y -n vibrant-env python=3
		#conda install -n vibrant-env -y -c bioconda prodigal hmmer
		#conda install -n vibrant-env -y -c ostrokach gzip
		#conda install -n vibrant-env -y -c conda-forge tar biopython matplotlib
		#conda install -n vibrant-env -y -c anaconda wget pandas seaborn numpy scikit-learn
		#conda activate vibrant-env
		#go_to_cloned_tools_folder
		#git clone https://github.com/AnantharamanLab/VIBRANT  
		#cp VIBRANT/databases/VIBRANT_setup.py $conda_path/envs/vibrant-env/bin
		#download-db.sh
		#pip install pickle-mixin
		#conda install -y -c bioconda vibrant==1.2.0
		#conda deactivate

		#INSTALLING stampede-clustergenomes
		#conda create -y -n stampede-clustergenomes-env 
		#conda install -n stampede-clustergenomes-env -y -c anaconda perl
		#conda install -n stampede-clustergenomes-env -y -c bioconda mummer
		#go_to_cloned_tools_folder
		#git clone https://bitbucket.org/MAVERICLab/stampede-clustergenomes.git
		#cd stampede-clustergenomes/bin
		#cp -r . $conda_path/envs/stampede-clustergenomes-env/bin
		#cd ../../..
		#conda activate stampede-clustergenomes-env
		#conda deactivate

		### to run stampede-clustergenomes is necessary to activate the env and run 
		### perl $your_path/mudoger_utils/cloned_tools/stampede-clustergenomes/bin/Stampede-clustergenomes.pl
		#cd ..
		
		#INSTALLING FASTA EXTRACTION ENV
		#conda create -y -n extract-env -c anaconda python=2.7.5
		#conda activate extract-env
		#conda deactivate

		#INSTALLING WISH
		#conda create -y -n wish-env
		#conda install -n wish-env -y -c conda-forge openmp
		#conda install -n wish-env -y -c anaconda make cmake
		#go_to_cloned_tools_folder
		#git clone https://github.com/soedinglab/WIsH.git
		#cd WIsH
		#conda activate wish-env
		#cmake .
		#make
		#conda deactivate wish-env
		#cp WIsH $conda_path/envs/wish-env/bin
		#cd ../..

		#INSTALLING CHECKV
		#conda create -y -n checkv-env -c conda-forge -c bioconda checkv

		#INSTALLING VCONTACT2
		#conda create -y -n vcontact2-env python=3 -c bioconda vcontact2
		#conda activate vcontact2-env
		#conda install -y -c bioconda vcontact2
		#conda install -y -c bioconda mcl blast diamond
		#conda install -y -c bioconda prodigal
		#conda update vcontact2
		## Install ClusterONE
		#go_to_cloned_tools_folder
		#wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar
		#cp cluster_one-1.0.jar $conda_path/envs/vcontact2-env/bin
		#cd .. #CHECK FILE LOCATION
		#conda deactivate

		break
	elif [ $choose = n -o $choose = N ]
	then
		echo "installation denied..."
		break
	else
		echo "command not found, please try again..."
	fi
done

echo "Do you want to install all the required tools for the Module 4 (Recovery of Eukaryotic Metagenome-Assembled Genomes)?[y/n]"
while :
do
	read choose
	if [ $choose = y -o $choose = Y ];
	then
		#INSTALLING EUKREP
		conda create -y -n eukrep-env -c bioconda eukrep

		#INSTALLING GENEMARKER-ES

		#INSTALLING MARKER2

		#INSTALLING BUSCO
		conda create -y -n busco-env -c bioconda busco

		#INSTALLING EUKCC
		conda create -y -n eukcc-env -c bioconda -c conda-forge eukcc

		break
	elif [ $choose = n -o $choose = N ]
	then
		echo "installation denied..."
		break
	else
		echo "command not found, please try again"
	fi
done

echo "Do you want to install all the required tools for the Module 5 (Relative Abundance)?[y/n]"
while :
do
	read choose
	if [ $choose = y -o $choose = Y ];
	then
		mkdir downloaded_tools && cd downloaded_tools
		conda create -y -n brat-env python=2.7
		conda activate brat-env

		#INSTALLING BRAT
		wget http://weaver.nlplab.org/~brat/releases/brat-v1.3_Crunchy_Frog.tar.gz
		tar xzf brat-v1.3_Crunchy_Frog.tar.gz
		cd brat-v1.3_Crunchy_Frog
		./install.sh -u
		python standalone.py

		conda deactivate
		cd ../..
		break
	elif [ $choose = n -o $choose = N ]
	then
		echo "installation denied..."
		break
	else
		echo "command not found, please try again"
	fi
done
