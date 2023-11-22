Installing MuDoGeR via conda can help the user to utilize only part of the workflow. However, it is recommended for those with a deeper understanding of how conda environments work, as manual adjustments may need to be made.

## Installation using Conda/Mamba
**1 - Install miniconda**

The MuDoGeR pipeline requires several tools that have multiple dependencies. More often than not, those tools require conflicting dependencies. To tackle this problem, MuDoGeR requires miniconda to be previously installed in your system. MuDoGeR will use miniconda to create multiple environments that are automatically orchestrated during the pipeline's functioning. Following you have a possible way of installing miniconda.

```console

#See documentation: https://docs.conda.io/en/latest/miniconda.html

$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

$ chmod +x Miniconda3-latest-Linux-x86_64.sh

$ ./Miniconda3-latest-Linux-x86_64.sh

$ export PATH=~/miniconda3/bin:$PATH

```

**2 - Install MuDoGeR**

Once you have miniconda installed and on your PATH, you can properly install MuDoGeR.
The MuDoGeR environment can be complex, and the installation script was designed to install and set up all necessary bioinformatics tools.

```console
#clone repository

$ git clone https://github.com/mdsufz/MuDoGeR.git

#Go to the MuDoGeR cloned repository folder
$ cd MuDoGeR

#Make sure you have conda ready and that you are in your base environment.
$ conda activate base
$ echo $CONDA_PREFIX

#You should see something like the following:
/path/to/miniconda3

#Install mamba in your base environment
$ conda install -c conda-forge mamba

#Run the installation script as follows
$ bash -i installation/install.sh

#Follow the instructions on the screen:
# Enter "y" if you want to install all modules, otherwise enter "n".
# If you entered "n",enter "y" for each of the modules you would like to install individually.

	The MuDoGeR's installation will begin..





	      (  )   (   )  )			
	       ) (   )  (  (			
	       ( )  (    ) )			
	       _____________			
	      <_____________> ___		
	      |             |/ _ \		
	      |               | | |		
	      |               |_| |		
	   ___|             |\___/		
	  /    \___________/    \		
	  \_____________________/		

	This might take a while. Time to grab a coffee...
```

**3 - Install necessary databases**

Several bioinformatics tools used within MuDoGeR require specific databases to work. We developed a database download and set up tool to make our lives easier. Make sure to run the database setup after MuDoGeR is installed.

You can also choose to install the databases used by each module individually. You can use the flag ```--dbs``` and choose the name of the module you want to set up the databases (all \[default], prokaryotes, viruses, eukaryotes).

Use this script if you want MuDoGeR to take care of everything. 

```console
#Make sure mudoger_env is activated. It should have been created when you ran 'bash -i installation/install.sh'
$ conda activate mudoger_env

#Go to MuDoGeR cloned directory
$ cd MuDoGeR

#Run database setup script
$ bash -i installation/database-setup.sh --dbs all -o /path/to/save/databases

#You can also check out the database-setup help information
$ bash -i installation/database-setup.sh --help

MuDoGeR database script v=1.0
Usage: bash -i database-setup.sh --dbs [module] -o output_folder_for_dbs

  --dbs all              		download and install the necessaries databases for all MuDoGeR modules [default]
  --dbs prokaryotes              	download and install the necessaries databases for prokaryotes module
  --dbs viruses              		download and install the necessaries databases for viruses module
  --dbs eukaryotes              	download and install the necessaries databases for eukaryotes module
  -o path/folder/to/save/dbs      	output folder where you want to save the downloaded databases
  --help | -h				show this help message
  --version | -v			show mudoger version


```

**3.1 - Update databases**

We plan to update the automatic database installation script at least once a year. The user can check the version of the database update script by using the ```--version``` flag. Finally, in case the user would like to manually update the databases, one should find the location of the databases used by MuDoGeR by running ```cat $CONDA_PREFIX/envs/mudoger_env/bin/database.sh```.

The user can then follow the instructions from the tool developer to update the desired database. The only requirement is that the root folder for the database uses the name of the tool as follows: ```buscodbs/  checkm/  checkv/  eukccdb/  gtdbtk/  vibrant/  wish/```.

**4 - Additional module 4 (eukaryotes) installation instructions**

Some tools used in module 4 (GENEMARK and MAKER2) require the user to provide information to the developers. Consequently, we could not implement an automatic installation and setup script. However, we created a tutorial to finish the module 4 setup.

The module 4 setup tutorial is found in [Module 4 setup](https://github.com/mdsufz/MuDoGeR/blob/master/installation/genemark_maker2_installation.md).

# MuDoGeR simplified usage - with conda environments installation

Currently, MuDoGeR v1.0 only works with paired-end ILLUMINA sequences. Future updates will add tools to work with long-read sequencing samples.
MuDoGeR was designed to work module by module, starting from pre-process (Module 1). Additional modularity will be added in future updates to allow the user to run specific parts of the pipeline. However, you can always use the tools independently by using the created conda environments by MuDoGeR. You can follow the instructions [here](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md#using-the-tools-independently).

**MuDoGeR** is an easy-to-use wrapper of several tools organized within modules. The individual modules can be called independently.

The pipeline requires, as input, a metadata table in tsv format containing the samples to be processed and the path to its raw sequence reads. The metadata file should have the sample name and the path to the forward reads file from the sample in one line, followed by the same sample name and the path to the reverse reads from the sample. An example metadata file is as follows:
```
#Show the content of the metadata.tsv file
$ cat metadata.tsv

EA_ERX4593008   /path/to/EA_ERX4593008/raw_reads_1.fastq
EA_ERX4593008   /path/to/EA_ERX4593008/raw_reads_2.fastq
EA_ERX4593009   /path/to/EA_ERX4593009/raw_reads_1.fastq
EA_ERX4593009   /path/to/EA_ERX4593009/raw_reads_2.fastq
EA_ERX4593010   /path/to/EA_ERX4593010/raw_reads_1.fastq
EA_ERX4593010   /path/to/EA_ERX4593010/raw_reads_2.fastq
EA_ERX4593011   /path/to/EA_ERX4593011/raw_reads_1.fastq
EA_ERX4593011   /path/to/EA_ERX4593011/raw_reads_2.fastq

```

## Please note that the forward sequencing reads file must end in "_1.fastq" and the reverse in "_2.fastq"! 

MuDoGeR is designed to run all multi-domain genome recovery pipelines entirely. In order for MuDoGeR to work automatically, from start to finish, we use a specific folder structure. Please read the [Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md) if you would like to manipulate MuDoGeR.

## MuDoGeR is installed via conda

Once MuDoGeR is installed, you can test it as follows:
```console
$ conda activate mudoger_env
$ mudoger --help


	███    ███ ██    ██ ██████   ██████   ██████  ███████ ██████  
	████  ████ ██    ██ ██   ██ ██    ██ ██       ██      ██   ██ 
	██ ████ ██ ██    ██ ██   ██ ██    ██ ██   ███ █████   ██████  
	██  ██  ██ ██    ██ ██   ██ ██    ██ ██    ██ ██      ██   ██ 
	██      ██  ██████  ██████   ██████   ██████  ███████ ██   ██ 
			Multi-Domain Genome Recovery
				Version 1.0.1



Mudoger v=1.0.1
Usage: mudoger --module module_name --meta metadata_table.tsv -o output_folder [module_options]

  --meta              		 Metadata table with your samples, as explained in the github documentation
  --module preprocess            Run all steps from module 1 (read_qc, kmer memory prediction and assembly)
  --module prokaryotes           Recovery of Prokaryotic Metagenome-Assembled Genomes
  --module viruses		 Recovery of Uncultivated Viral Genomes
  --module eukaryotes		 Recovery of Eukaryotic Metagenome-Assembled Bins
  --module abundance_tables	 pMAGs/UViGs/eMABs coverage and abundance calculation
          	type             can be --reduced (default) , --complete or --genes
          	mapping_type	 can be --absolute-values (default), --coverage, and --relative-abundance
  --help | -h		         show this help message
  --version | -v	         show mudoger version

```

A simplified use of MuDoGeR can be done as follows:

```console

$ mudoger --module preprocess --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20 -m 100

$ mudoger --module prokaryotes --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20

$ mudoger --module viruses --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20

$ mudoger --module eukaryotes --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20

$ mudoger --module abundance_tables --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20 --reduced --coverage --relative-abundance

```
