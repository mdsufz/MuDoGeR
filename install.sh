#!/bin/bash

#making a directory for mudoger and change directory to it
mkdir mudoger
cd mudoger

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
		mkdir cloned_tools
		cd cloned_tools
		git clone https://github.com/bxlab/metaWRAP.git
		#->configure the config-metawrap file to pointo to the database
		conda create -y -n metawrap-env python=2.7
		cd metaWRAP/bin
		cp -r . $conda_path/envs/metawrap-env/bin
		cd ../../..
		conda activate metawrap-env
		conda config --add channels defaults
		conda config --add channels conda-forge
		conda config --add channels bioconda
		conda config --add channels ursky
		conda install -y --only-deps -c ursky metawrap-mg=1.3.2
		metawrap
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
		#INSTALLING GTDB-TK
		conda create -y -n gtdbtk-env
		conda install -n gtdbtk-env -y -c bioconda gtdbtk

		#INSTALLING UBIN
		cd cloned_tools || mkdir cloned_tools && cd cloned_tools
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
		#INSTALLING VIRSORTER
		conda create --name virsorter-env -y -c bioconda mcl=14.137 muscle blast perl-bioperl perl-file-which hmmer=3.1b2 perl-parallel-forkmanager perl-list-moreutils diamond=0.9.14
		cd cloned_tools || mkdir cloned_tools && cd cloned_tools
		git clone https://github.com/simroux/VirSorter.git
		cd VirSorter/Scripts
		conda install -n virsorter-env -y -c anaconda make
		conda activate virsorter-env
		make clean
		make
		conda deactivate
		cp -r . $conda_path/envs/virsorter-env/bin
		cd ..
		cp wrapper_phage_contigs_sorter_iPlant.pl $conda_path/envs/virsorter-env/bin
		conda install -n virsorter-env -y -c bioconda metagene_annotator
		cd ../..

		#INSTALLING VIRFINDER
		conda create -y -n virfinder-env
		conda install -n virfinder-env -y -c bioconda r-virfinder

		#INSTALLING VIBRANT
		conda create -y -n vibrant-env python=3
		conda install -n vibrant-env -y -c bioconda prodigal hmmer
		conda install -n vibrant-env -y -c ostrokach gzip
		conda install -n vibrant-env -y -c conda-forge tar biopython matplotlib
		conda install -n vibrant-env -y -c anaconda wget pandas seaborn numpy scikit-learn
		conda activate vibrant-env
		pip install pickle-mixin
		conda install -y -c bioconda vibrant==1.2.0
		conda deactivate

		#INSTALLING stampede-clustergenomes
		conda create -y -n stampede-clustergenomes-env 
		conda install -n stampede-clustergenomes-env -y -c anaconda perl
		conda install -n stampede-clustergenomes-env -y -c bioconda mummer
		cd cloned_tools || mkdir cloned_tools && cd cloned_tools
		git clone https://bitbucket.org/MAVERICLab/stampede-clustergenomes.git
		#to run stampede-clustergenomes is necessary to activate the env and run perl $your_path/mudoger/cloned_tools/stampede-clustergenomes/bin/Stampede-clustergenomes.pl
		cd ..

		#INSTALLING WISH
		conda create -y -n wish-env
		conda install -n wish-env -y -c conda-forge openmp
		conda install -n wish-env -y -c anaconda make cmake
		cd cloned_tools || mkdir cloned_tools && cd cloned_tools
		git clone https://github.com/soedinglab/WIsH.git
		cd WIsH
		conda activate wish-env
		cmake .
		make
		conda deactivate wish-env
		cp WIsH $conda_path/envs/wish-env/bin
		cd ../..

		#INSTALLING CHECKV
		conda create -y -n checkv-env -c conda-forge -c bioconda checkv

		#INSTALLING VCONTACT2
		conda create -y -n vcontact2-env -c bioconda vcontact2
		conda activate vcontact2-env
		conda update vcontact2
		conda deactivate

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