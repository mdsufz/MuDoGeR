 # Multi-Domain Genome Recovery v1.0.1 (MuDoGeR v1.0.1)
 
 
![ScreenShot](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/figure_1_overview.jpg)


The **Multi-Domain Genome Recovery v1.0 (MuDoGeR v1.0)** framework ([**Figure 1**](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/figure_1_overview.jpg)) is a tool developed to help users to recover Metagenome-Assembled Genomes (MAGs as defined by [Parks et al. (2018)](https://www.nature.com/articles/s41564-017-0012-7)) and Uncultivated Viral Genomes (UViGs as defined by  [Roux (2019)](https://www.nature.com/articles/nbt.4306)) from whole-genome sequence (WGS) samples simultaneously. The **MuDoGeR v1.0** framework acts as a wrapper for several tools. It was designed to be an easy-to-use tool that outputs ready-to-use comprehensive files.

You should be able to run 1 simple command for each module. Therefore, you only need 5 commands to completely run MuDoGeR. After a successful run of MuDoGeR, you should have the outputs summarized in [Figure 2](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/mudoger_outputs.jpg). Please find a comprehensive description of the main outputs in the [understand outputs](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md) file.

![ScreenShot](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/mudoger_outputs.jpg)

### Using the tools individually

In addition, MuDoGeR also sets up individual working conda environments for each of the integrated tools. Consequently, if the user wants to customize the use of any tool, you can use MuDoGeR to configure your machine and follow the instructions [here](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md#using-the-tools-independently) to activate the relevant environments.

## Reading this GitHub

This Github should help you install and run the complete MuDoGeR pipeline, as well as understand all its outputs. Consequently, we suggest the following reading strategy: 

* First, read the [MuDoGeR overview](https://github.com/mdsufz/MuDoGeR#mudoger-overview) and define which modules you are interested in using.
* Secondly, read the [System requirements](https://github.com/mdsufz/MuDoGeR#system-requirements) and make sure you have the resources for the modules you want to use.
* Then, read the [Installation](https://github.com/mdsufz/MuDoGeR#installation) and follow its steps.
* Read the [overview descriptions](https://github.com/mdsufz/MuDoGeR#modules-overview) of the MuDoGeR modules you intend to use.
* If you want a quick run, read the [MuDoGeR simplified usage](https://github.com/mdsufz/MuDoGeR#mudoger-simplified-usage)
* Read the [understanding main outputs](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md) file. To understand and find relevant outputs created by MuDoGeR.
* Read the [MuDoGeR Manual](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md) for a more detailed description of the used modules and their output files

## MuDoGeR Overview

**MuDoGeR** was designed to be an easy-to-use tool that outputs ready-to-use comprehensive files. You should be able to run 1 simple command for each of the following modules.

The **MuDoGeR** starts with Module 1: **Pre-Processing**, which covers the Raw Read Quality Control, Resources calculation, and Assembly. The assembled libraries should be used in all the other modules.

After the data pre-processing, **MuDoGeR** is divided into 3 different branches:
Module 2: **Recovery of Prokaryotic Metagenome-Assembled Genomes (pMAGs)**
Module 3: **Recovery of Uncultivated Viral Genomes (UViGs)**
Module 4: **Recovery of Eukaryotic Metagenome-Assembled Bins (eMABs)**

Furthermore, in **Module 5**: **Relative Abundance**, users can automatically calculate the coverage and relative abundance tables from the recovered pMAGs/UViGs/eMABs. The users can also calculate the coverage and relative abundance tables from the prokaryotic genes annotated in assembled libraries. 

 
* Please find a comprehensive description of the main outputs in the [understand outputs](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md) file.
* Instructions for using the **MuDoGeR** can be found in the following hyperlink: [Manual MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md).
* Information about the system requirements of the **MuDoGeR** can be found in the following hyperlink: [System requirements](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#system-requirements).
* Detailed instructions for the installation of the **MuDoGeR** tools can be found in the following hyperlink: [Installation](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#installation).
* The simplified usage of the **MuDoGeR** can be found in the following hyperlink: [MuDoGeR simplified usage](https://github.com/mdsufz/MuDoGeR#mudoger-simplified-usage).
* To use the individual working conda environments created by MuDoGeR for each of the used tools, go [here](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md#using-the-tools-independently).


# System requirements

**MuDoGeR makes it easy to install and run a group of complex software and dependencies.** More than 20 bioinformatics tools and several complex software dependencies are present in MuDoGeR. Hopefully, you won't have to worry too much about it. We designed an installation script that should take care of every dependency for you. You can find the **MuDoGeR** installation tutorial [here](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#installation).

Keep in mind that the MuDoGeR pipeline requires some computer power, and you probably won't be able to run it on a laptop. The complete software installation requires approximately 170 GB, but **MAKER2**, from **Module 4**, uses 99 GB of that space since it requires the database to be installed in a specific manner. See [Module 4 setup](https://github.com/mdsufz/MuDoGeR/blob/master/installation/genemark_maker2_installation.md). The complete database requirements, considering all tools, is around 439.9 GB. However, you don't need to install all of MuDoGeR's Modules to use it.

MuDoGeR is designed to support only Linux x64 systems. As for the resource requirements, the MuDoGeR framework uses software that requires a large amount of RAM (e.g. **GDTB-Tk**, **MetaWRAP** ). Specific resource requirements vary depending on your data and its sequencing depth. We recommend the user provide at least 180 GB of RAM. 
Therefore, for the assembly process, **MuDoGeR** attempts to calculate the amount of memory necessary for **metaSPades** (on step **1.b**). The user should be aware that samples with higher expected diversity require a higher amount of memory.

Consequently, we suggest you install and run **MuDoGeR** using your available high-performance computer or in cloud services such as [AWS](https://aws.amazon.com/ec2/?nc2=h_ql_prod_fs_ec2), [Google cloud](https://cloud.google.com/compute#section-8), or, for researchers in Germany, the [de.NBI](https://www.denbi.de/)

The software dependencies used during the pipeline are described here: [Dependencies description](https://github.com/mdsufz/MuDoGeR/blob/master/dependencies_description.md).

# Installation

## Installation using Singularity (now called Apptainer) - Recommended

**0 - Install Singularity**
Most HPC administrators are already making Singularity available for its users. You could check if that is your case and skip this step.
Otherwise, please follow the instructions on the official [Singularity installation guide](https://apptainer.org/admin-docs/master/installation.html)

**1 - Download MuDoGeR ready-to-use container**

Once you have Singularity installed, you can simply download the MuDoGeR container. Remember that the container's usage is slightly different. Please refer to the Singularity container [usage](https://github.com/mdsufz/MuDoGeR#mudoger-singularity-usage-notes)

Click [here](https://e.pcloud.link/publink/show?code=XZ57q2ZksirrSUwBT0wEuorJjhuih93vQnV) to be redirected to the download page.

Following, you can click on Direct Download or right-click it and "copy link". Once with the copied link, you can use ```wget``` on your platform.


**2 - Database installation**

The MuDoGeR required databases can vary depending on which module you plan to use. Naturally, the databases can require significant storage and are not included in the MuDoGeR container.
The user can then follow the instructions from the tool developer to install and update the desired database.
The only requirement is that all the databases use the same base folder and are installed using the name of the tool as follows: ```buscodbs/  checkm/  checkv/  eukccdb/  gtdbtk/  vibrant/  wish/```.
Therefore, your database installation folder should look like this:

```
mudoger_dbs/
├── buscodbs
├── checkm
├── checkv
├── eukccdb
├── gtdbtk
├── vibrant
└── wish
```

**3 - Configure Genemark License if you will use Module 4**

1. ACCESS GENEMARK WEBPAGE

  http://opal.biology.gatech.edu/GeneMark/license_download.cgi


2. SELECT OPTIONS 

      **GeneMark-ES/ET/EP ver \*_lic and LINUX 64**


3. FILL IN THE CREDENTIALS WITH ***YOUR NAME, E-MAIL, INSTITUTION, ETC... ***


4. CLICK ON ***'I agree the terms of this license agreement'***


5. DOWNLOAD THE 64_bit key files provided

      It should look something like the following:

        ```console
        $ wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_HZzc0/gm_key_64.gz
        ```

      You should have the following file:  **gm_key_64.gz**


6. DECOMPRESS THE KEY FILE

```console
$ gunzip gm_key_64.gz
```

7. COPY AND RENAME KEY FILE TO A FOLDER

The folder you will move the renamed key file will be used as your Home during the execution of Module 4 in singularity. Please see [here](https://github.com/mdsufz/MuDoGeR#mudoger-singularity-usage-notes) to run Module 4 using the singularity container. 

```console
$ cp gm_key_64 /path/to/folder/.gm_key
```

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

The module 4 setup tutorial is found in [Module 4 setup](https://github.com/mdsufz/MuDoGeR/blob/master/installation/genemark_maker2_installation.md)

## Modules Overview
### Module 1: Pre-Processing 

![Screenshot](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/module_1_ovrw.jpg) 

The steps of Module 1 are shown in [Figure 3](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/module_1_ovrw.jpg). A detailed description of its execution and outputs are found here: [Pre-Processing description](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-1-pre-processing).
 
 When you use MuDoGeR Module 1, it will perform the following tasks:

* **1.a**: Raw Read Quality Control.
* **1.b**: Calculation of memory requirements for the assembly process.
    * **(1.b.1)** The k-mer (33-mer and 55-mer) of the quality-controlled reads produced in **1.a** is calculated.  
    * **(1.b.2)** The calculated k-mer is used in a trained machine learning model to estimate the amount of memory that **metaSPades** uses to assemble the reads.
* **1.c**: Assembly of the quality-controlled reads.

### Module 2: Recovery of Prokaryotic Metagenome-Assembled Genomes (pMAGs)

![screenshot](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/module_2_ovrw.jpg)

Module 2 workflow is shown in [Figure 4](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/module_2_ovrw.jpg). A detailed description of its execution and outputs are found here: [Pipeline for recovery of Prokaryotic Metagenome-Assembled Genomes](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-2-recovery-of-prokaryotic-metagenome-assembled-genomes-pmags).

When you use MuDoGeR Module 2, it will perform the following tasks:
 
* **2.a**: Binning and bin refinement of the Prokaryotic bins.
     * **(2.a.1)** Binning with [**Metabat2**](https://peerj.com/articles/7359/), [**Maxbin2**](https://academic.oup.com/bioinformatics/article/32/4/605/1744462), and [**CONCOCT**](https://www.nature.com/articles/nmeth.3103).
     * **(2.a.2)** Bacterial bins refinement and archaea bins refinement.
     * **(2.a.3)** Dereplication of the recovered prokaryotic bins.
* **2.b**: Taxonomic classification, quality estimation, and gene annotation.
     * **(2.b.1)** Taxonomic classification of the prokaryotic bins produced in **(2.a.3)** using [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182).
     * **(2.b.2)** Generation of quality matrix of the prokaryotic bins produced in **(2.a.3)** using [**CheckM**](https://genome.cshlp.org/content/25/7/1043). 
     * **(2.b.3)** Prokaryotic MAGs gene annotation with [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517).
* **2.c**: Sequence metrics calculation and selection of Prokaryotic MAGs.
     * **(2.c.1)** Sequence metric calculation from the selected MAGs.
     * **(2.c.2)** Selection of Prokaryotic MAGs
    
    
### Module 3: Recovery of Uncultivated Viral Genomes (UViGs)

![screenshot](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/module_3_ovrw.jpg) 

The steps of Module 3  are shown in [Figure 5](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/module_3_ovrw.jpg). A detailed description of its execution and outputs are found here: [Pipelines for viral genomes recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-3-recovery-of-uncultivated-viral-genomes-uvigs).

When you use MuDoGeR Module 3, it will perform the following tasks: 
 
* **3.a**: Recovery of Putative Viral Contigs
    * **(3.a.1)** Recovery of putative viral contigs using [**VirSorter2**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y), [**VirFinder**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0283-5), and [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0).
    * **(3.a.2)** Filtering of the putative viral contigs.
    * **(3.a.3)** Dereplication of the putative viral contigs.
 * **3.b**: Taxonomic and Quality estimation of potential viral contigs
    * **(3.b.1)** Taxonomic classification from the dereplicated putative viral contigs with [**Vcontact2**](https://www.nature.com/articles/s41587-019-0100-8).
    * **(3.b.2)** Checking the quality of the dereplicated contigs with [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7).
* **3.c**: Viral-Host pair estimation using [**WIsH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377). This step is only done automatically if you generate the prokaryotic MAGs using MuDoGeR as well. 
* **3.d**: Selection of UViGs
    * **(3.d.1)** Selection of all viruses that yielded taxonomy when using vContact2 plus those larger than 15 Kb.
    * **(3.d.2)** Selection based on the quality determined by [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7).

### Module 4: Recovery of Eukaryotic Metagenome-Assembled Bins (eMABs)

![screenshot](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/module_4_ovrw.jpg)

The steps of Module 4  are shown in [Figure 6](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/module_4_ovrw.jpg). A detailed description of its execution and outputs are found here:  [Pipelines for eukaryotic bins recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-4-recovery-of-eukaryotic-metagenome-assembled-bins-emabs).

When you use MuDoGeR Module 4, it will perform the following tasks:

* **4.a**: Recovery and binning of Eukaryotic assemblies.
    * **(4.a.1)** Classification of Eukaryotic assemblies and removal of prokaryotic assemblies with [**EukRep**](https://genome.cshlp.org/content/28/4/569).
    * **(4.a.2)** Use of [**CONCOCT**](https://www.nature.com/articles/nmeth.3103) for binning the Eukaryotic assemblies. 
    * **(4.a.3)** Filtering the Eukaryotic bins, produced from [**CONCOCT**](https://www.nature.com/articles/nmeth.3103), by size. Bins with size < 1.5 Mb are removed.
* **4.b**: Completeness and contamination estimation and annotation of Eukaryotic bins
    * **(4.b.1)** In the filtered bins produced in **4.a**, genes are predicted using [**GeneMark**](https://academic.oup.com/nar/article/29/12/2607/1034721).
    * **(4.b.2)** Completeness and contamination estimation of the Eukaryotic filtered bins produced in **4.a** using [**EukCC**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02155-4).
    * **(4.b.3)** [**MAKER2**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491) annotates the predicted genes produced by [**GeneMark**](https://academic.oup.com/nar/article/29/12/2607/1034721). 
    * **(4.b.4)** [**BUSCO**](https://academic.oup.com/bioinformatics/article/31/19/3210/211866) is applied to the annotated genes from [**MAKER2**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491), for detection of single-copy orthologous genes (SCGs) and estimation of completeness of Eukaryotic contigs.


### Module 5 Relative abundance

![screenshot](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/module_5_ovrw.jpg)

The steps of Module 5  are shown in [Figure 7](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/module_5_ovrw.jpg). A detailed description of its execution and outputs are found here: [Pipelines for abundance calculation](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-5-relative-abundance).
Essentially, module 5 maps the quality-controlled reads of your sample on the recovered pMAGs/UViGs/eMAB or annotated prokaryotic genes. We designed three possible mapping types to calculate abundance: **reduced**, **complete**, or **genes**. A detailed description of their differences can be found [here](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-5-relative-abundance)

The steps of Module 5 can be summarized as follows. If you select **complete** or **genes**, the pipeline will run steps **5.a** and **5.b**. If you select **genes**, the pipeline will run **5.c**: 

* **5.a**: Select representative pMAGs from each created OTU
    * **(5.a.1)** Copy recovered pMAGs from all samples within the provided metadata table.
    * **(5.a.2)** Group recovered pMAGs of all samples within the provided metadata table using **gOTUpick**.
    * **(5.a.3)** Select the highest quality pMAG within the **gOTUpick** groups as the group's representative MAG.
* **5.b**: pMAGs/UViGs/eMABs mapping and abundance calculation
    * **(5.b.1)** Copy representative pMAGs from step **5.1** and the selected UViGs and eMABs recovered in module 3 and module 4, respectively.
    * **(5.b.2)** Index selected pMAGs/UViGs/eMABs
    * **(5.b.3)** If ```--coverage``` is selected, calculate pMAGs/UViGs/eMAB size, and the average read length of all samples to be mapped. If ```--relative-abundance``` calculate the total number of reads from all samples. This information is used further in the pipeline.
    * **(5.b.4)** If ```--reduced``` is selected, maps reads from the samples where the pMAGs/UViGs/eMABs were found on the pMAGs/UViGs/eMABs. If ```--complete``` is selected, map reads from all samples on the pMAGs/UViGs/eMABs.
    * **(5.b.5)** Calculate the absolute number of hits, relative abundance, and coverage tables, if the respective flag is selected.
* **5.c**: Genes relative abundance calculation from the samples assembly. Currently working on prokaryotic genes
    * **(5.c.1)** Index assemblies from given samples.
    * **(5.c.2)** Map sample reads on the respective assembly.
    * **(5.c.3)** Annotate genes on the assembly with [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517).
    * **(5.c.4)** Convert the .gff file from [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) gene annotation into .gtf.
    * **(5.c.5)** Count mapped reads on each gene.
    * **(5.c.6)** Calculate the average read length of all mapped samples.
    * **(5.c.7)** Calculate gene lenght from [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) .gtf file.
    * **(5.c.8)** Calculate genes' absolute number of hits, relative abundance, coverage, and TPM tables for each sample.


# MuDoGeR simplified usage

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

MuDoGeR is designed to run all multi-domain genome recovery pipelines entirely. In order for MuDoGeR to work automatically, from start to finish, we use a specific folder structure. Please, read the [Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md) if you would like to manipulate MuDoGeR. 

## MuDoGeR Singularity Usage Notes

When using the MuDoGeR singularity container, you have all the complex dependencies and software environments from MuDoGeR already configured.
To keep things simple, we suggest you keep your metadata.tsv file in the same folder where you have your samples.
To run the container, you should use the ```--bind``` flag from singularity to mount your data to the container. For instance, if your input data and metadata.tsv is located at ```/path/to/input```, your databases are in ```/path/to/mudoger_dbs```, and you want to save the results in ```/path/to/output```, your command call will look as follows:

**Do not change the folders ```/tools/data_input, /tools/mudoger_output_in_container, /tools/dbs, and /mudoger_home ```**, as those are the folder used inside the container.

One additional point of attention. The input data is mounted in the container at ```/tools/data_input/```. Therefore, if your sample sequence files are in ```/path/to/input/sampleID/sampleID_1.fastq and /path/to/input/sampleID/sampleID_2.fastq``` your metadata.tsv file should look like:

```
#Show the content of the metadata.tsv file
$ cat metadata.tsv

sampleID   /tools/data_input/sampleID/sampleID_1.fastq
sampleID   /tools/data_input/sampleID/sampleID_2.fastq
```
Usage:
```console

#singularity exec --bind /path/to/input:/tools/data_input,/path/to/output:/tools/mudoger_output_in_container,/path/to/mudoger_dbs:/tools/dbs mudogerV1.sif mudoger [COMMAND]

singularity exec --bind /path/to/input:/tools/data_input,/path/to/output:/tools/mudoger_output_in_container,/path/to/mudoger_dbs:/tools/dbs mudogerV1.sif mudoger --module preprocess --meta /tools/data_input/metadata.tsv -o /tools/mudoger_output_in_container -t 25 -m 100

singularity exec --bind /path/to/input:/tools/data_input,/path/to/output:/tools/mudoger_output_in_container,/path/to/mudoger_dbs:/tools/dbs mudogerV1.sif mudoger --module prokaryotes --meta /tools/data_input/metadata.tsv -o /tools/mudoger_output_in_container -t 25

singularity exec --bind /path/to/input:/tools/data_input,/path/to/output:/tools/mudoger_output_in_container,/path/to/mudoger_dbs:/tools/dbs mudogerV1.sif mudoger --module viruses --meta /tools/data_input/metadata.tsv -o /tools/mudoger_output_in_container -t 25

singularity exec --bind /path/to/input:/tools/data_input,/path/to/output:/tools/mudoger_output_in_container,/path/to/mudoger_dbs:/tools/dbs mudogerV1.sif mudoger --module abundance_tables --meta /tools/data_input/metadata.tsv -o /tools/mudoger_output_in_container -t 25 --reduced --coverage --relative-abundance

```

Module 4 (Eukaryotes recovery) has one particularity. The [**GeneMark**](https://academic.oup.com/nar/article/29/12/2607/1034721) requires each user to agree to a license and place it in their home folder. You can obtain this license following the instructions [here](https://github.com/mdsufz/MuDoGeR#installation-using-singularity-now-called-apptainer---recommended). Once the license is configured, you have to ```--bind``` the folder containing it to the singularity container using the ```--home``` flag. For instance, if you saved the Genemark license key in ```/path/to/tmp_home``` your module 4 singularity command will be:

```console

singularity exec --bind /path/to/input:/tools/data_input,/path/to/output:/tools/mudoger_output_in_container,/path/to/mudoger_dbs:/tools/dbs,/path/to/tmp_home:/mudoger_home --home /mudoger_home mudogerV1.sif mudoger --module eukaryotes --meta /tools/data_input/metadata.tsv -o /tools/mudoger_output_in_container/ -t 25

```



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

The result two-level folder structure after a successful run of all MuDoGeR is as follows:
```console
.
├── sample_1
│   ├── assembly
│   ├── eukaryotes
│   ├── khmer
│   ├── prokaryotes
│   ├── qc
│   └── viruses
├── sample_2
│   ├── assembly
│   ├── eukaryotes
│   ├── khmer
│   ├── prokaryotes
│   ├── qc
│   └── viruses
└── mapping_results
    ├── assembly_gene_map
    ├── euk_mabs_mapping
    ├── gOTUpick_results
    ├── merged_reads
    ├── pmags_otu_mapping
    └── uvigs_mapping
```

A more detailed tutorial for the MuDoGeR can be found in [Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md).

# MuDoGeR as Wrapper and its Critical Use

MuDoGeR is a wrapper designed to streamline the genome assembly process from metagenome samples across multiple domains. While MuDoGeR accelerates metagenomics analysis, it's crucial to understand the inherent limitations of any metagenomic approach. The tools and parameters integrated into MuDoGeR are based on benchmark studies, but users should understand their dataset and tools' limitations to adapt the workflow.

* MuDoGeR generates progress reports, intermediate files, and error logs at various stages of the assembly process, which are detailed in the [MuDoGeR Manual](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md). Users should regularly check these reports to ensure the tool works as expected for their dataset.

* MuDoGeR aims to help provide a holistic view of all three domains simultaneously, reducing cross-domain recovery bias by initiating all genome recovery from the same assembly. However, the genome recovery approaches from different domains are at different stages of technological progress and complexity of analysis.
* Users should be aware of potential recovery bias from a particular dataset.
* For experienced users, MuDoGeR allows the activation of each tool separately, enabling users to adapt the process to their specific needs. follow the instructions [here](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md#using-the-tools-independently).
* Users are strongly recommended to consult and check the direct links from the tools used within the wrapper for a deeper understanding of the underlying processes and optimization of parameters for their specific needs. The software used during the pipeline is detailed and described here: [Dependencies description](https://github.com/mdsufz/MuDoGeR/blob/master/dependencies_description.md).

# Citing

A preprint of the manuscript can be found in [Biorxiv](https://www.biorxiv.org/content/10.1101/2022.06.21.496983v3.full)

# Acknowledgements

	

 


	

