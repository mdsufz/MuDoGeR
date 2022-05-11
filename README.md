 # Multi-Domain Genome Recovery (MuDoGeR)
 
 
![ScreenShot](https://github.com/JotaKas/MuDoGeR/blob/master/flowcharts/fig1_20.5.21.png)



The **Multi-Domain Genome Recovery (MuDoGeR)** framework (**Figure 1**) is a tool developed to help users to recover Metagenome-Assembled Genomes (MAGs as defined by Parks et al. (2018)) and Uncultivated Viral Genomes (UViGs as defined by  Roux (2019)) from whole-genome sequence (WGS) samples simultaneously. The **MuDoGeR** framework act as a wrapper of several tools. It was designed to be an easy-to-use tool that outputs ready-to-use comprehensive files.


The **MuDoGeR** starts with Module 1: **Pre-Processing**, which covers: **1.a** **Raw Read Quality Control** and **1.b** **Resources calculation** and **1.c** **Assembly**. The assembled sequences should be used in all the other modules.

After pre-processing of the data, **MuDoGeR** is divided into 3 different branches:
Module 2: **Recovery of Prokaryotic Metagenome-Assembled Genomes**
Module 3: **Recovery of Uncultivated Viral Genomes**
Module 4: **Recovery of Eukaryotic Metagenome-Assembled Genomes**

Furthermore, in **Module 5**: **Relative Abundance**, there are the selection of the MAGs from the group of processed WGS samples, and the calculation of the coverage and relative abundance table of the selected MAGs and UViGs within the sample's group. References of the used tools can be found at the end of the page.

## MuDoGeR Overview

* Instructions for using the **MuDoGeR** can be found in the following hyperlink: [Manual MuDoGeR](https://github.com/JotaKas/MuDoGeR/blob/master/Manual_MuDoGeR.md).
* Detailed description of the **MuDoGeR** steps can be found in the following hyperlink: [Module description](https://github.com/JotaKas/MuDoGeR/blob/master/module_description.md).
* Information about the system requirements of the **MuDoGeR** can be found in the following hyperlink: [System requirements](https://github.com/JotaKas/MuDoGeR/blob/master/README.md#system-requirements).
* Detailed instructions for the installation of the **MuDoGeR** tools can be found in the following hyperlink: [Installation](https://github.com/JotaKas/MuDoGeR/blob/master/README.md#installation).
* The simplified usage of the **MuDoGeR** can be found in the following hyperlink: [Simplified usage of the MuDoGeR](https://github.com/JotaKas/MuDoGeR/blob/master/README.md#simplified-usage-of-the-mudoger).
### Module 1: Pre-Processing 

![Screenshot](https://github.com/JotaKas/MuDoGeR/blob/master/flowcharts/Module%201.PNG) 

The steps of Module 1 are shown in Figure 2. A detailed description of its execution and outputs are found here: [Pre-Processing description](https://github.com/JotaKas/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-1-pre-processing).
 
 The steps of Module 1 can be summarized as follows:

* **1.a**: Raw Read Quality Control.
* **1.b**: Calculation of memory requirements for the assembly process.
    * **(1.b.1)** The k-mer (33-mer and 55-mer) of the quality-controled reads produced in **1.a** is calculated.  
    * **(1.b.2)** The calculated k-mer is used within and trained machine learning model to estimate the amount of memory that **metaSPades** use to assemble the reads.
* **1.c**: Assembly of the quality-controlled reads.

### Module 2: Recovery of Prokaryotic Metagenome-Assembled Genomes (MAGs)

![screenshot](https://github.com/JotaKas/MuDoGeR/blob/master/flowcharts/MODULE_2_21.5.21.png)

Module 2 workflow is shown in Figure 3. A detailed description of its execution and outputs are found here: [Pipeline for recovery of Prokaryotic Metagenome-Assembled Genomes](https://github.com/JotaKas/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-2-recovery-of-prokaryotic-metagenome-assembled-genomes).

 The steps of Module 2 can be summarized as following:
 
* **2.a**: Binning and bin refinement of the Prokaryotic bins.
     * **(2.a.1)** Binning with **MaxBin2**, **metaBAT2**, and **CONCOCT**.
     * **(2.a.2)** Bacterial bins refinement and archea bins refinement.
     * **(2.a.3)** Dereplication of the recovered prokaryotic bins.
* **2.b**: Taxonomic classification, quality estimation, and gene annotation. selection of Prokaryotic MAGs.
     * **(2.b.1)** Taxonomic classification of the prokaryotic bins produced in **(2.a.3)** using **GTDB-Tk**.
     * **(2.b.2)** Generation of quality matrix of the prokaryotic bins produced in **(2.a.3)** using **CheckM**. 
     * **(2.b.3)** Prokaryotic MAGs gene annotation with **PROKKA**.
* **2.c**: Sequence metrics calculation and selection of Prokaryotic MAGs.
     * **(2.c.1)** Sequence metric calculation from the selected MAGs.
     * **(2.c.2)** Prokaryotic MAGs gene annotation with **PROKKA**.
    
    
### Module 3: Recovery of Uncultivated Viral Genomes (Uvigs)

![](https://github.com/JotaKas/MuDoGeR/blob/c65d851f6439cc4eb8c18672afb5cde9a9165f40/flowcharts/Module%203.png) 

The steps of Module 3  are shown in Figure 4. A detailed description of its execution and outputs are found here: [Pipelines for viral genomes recovery](https://github.com/JotaKas/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-3-recovery-of-uncultivated-viral-genomes).

 The steps of Module 3 can be summarized as follows: 
 
* **3.a**: Recovery of Uncultivated Viral Genomes
    * **(3.a.1)** Recovery of Uvigs using **VirFinder**, **VirSorter** and **VIBRANT**.
    * **(3.a.2)** Filtering of the Uvigs.
    * **(3.a.3)** Dereplication of the Uvigs.
 * **3.b**: Taxonomic and Quality estimation of Uvigs
    * **(3.b.1)** Taxonomic classification from the dereplicated Uvigs with **vContact2**.
    * **(3.b.2)** Checking the quality of the dereplicated contigs with **CheckV**.
* **3.c**: Host identification of the dereplicated Uvigs using **WIsH**. This step is only done automatically if you generate the prokaryotic MAGs using MuDoGeR as well. 
* **3.d**: Selection of Uvigs
    * **(3.d.1)** Selection of all viruses that yielded taxonomy when using vContact2 plus those larger than 15 Kb.
    * **(3.d.2)** Selection based on the quality determined by **CheckV**


### Module 4: Recovery of Eukaryotic Metagenome-Assembled Genomes (eMAGs)

[](https://github.com/JotaKas/MuDoGeR/blob/master/Module%204.png)

The steps of Module 4  are shown in Figure 5. A detailed description of its execution and outputs are found here:  [Pipelines for eukaryotic genomes recovery](https://github.com/JotaKas/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-4-recovery-of-eukaryotic-metagenome-assembled-genomes).

 The steps of Module 4 can be summarized as follows:
* **4.a**: Recovery and binning of Eukaryotic assemblies.
    * **(4.a.1)** Classification of Eukaryotic assemblies and removal of prokaryotic assemblies with **EukRep**.
    * **(4.a.2)** Use of **CONCOCT** for binning the Eukaryotic assemblies.  
    * **(4.a.3)** Filtering the Eukaryotic bins, produced from **CONCOCT**, by size. Bins with size < 1.5 Mb are removed.
* **4.b**: Completeness and contamination estimation and annotation of Eukaryotic bins
    * **(4.b.1)** In the filtered bins produced in **4.a**, genes are predicted using **GeneMark**.
    * **(4.b.2)** Completeness and contamination estimation of the Eukaryotic filtered bins produced in **4.a** using **EukCC**.
    * **(4.b.3)** **MAKER2** annotates the predicted genes produced by **GeneMark**. 
    * **(4.b.4)** **BUSCO** is applied to the annotated genes from **MAKER2**, for detection of single-copy orthologous genes (SCGs) and estimation of completeness of Eukaryotic contigs.


### Module 5 Relative abundance
[](https://github.com/JotaKas/MuDoGeR/blob/master/Module%205.PNG)

The steps of Module 5  are shown in Figure 6. A detailed description of its execution and outputs are found here: [Pipelines for abundance calculation](https://github.com/JotaKas/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-5-relative-abundance).
Essentially, module 5 maps the quality-controlled reads of your sample on the recovered MAGs or annotated prokaryotic genes. We designed three possible mapping types to calculate abundance: **reduced**, **complete**, or **genes**. A detailed description of their differences can be found [here](https://github.com/JotaKas/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-5-relative-abundance)

The steps of Module 5 can be summarized as follows. Depending on your choice, your samples will be mapped following steps on **5.a**, **5.b** , or **5.c**: 

* **5.a**: Complete relative abundance calculation
    * **(5.a.1)** Select MAGs from all samples within the provided metadata table.
    * **(5.a.2)** Calculate MAGs size, the number of quality-controlled reads and the average read length of all samples to be mapped.
    * **(5.a.3)** Index MAGs.
    * **(5.a.4)** Map reads from samples on their respective recovered MAGs.
    * **(5.a.5)** Calculate MAGs absolute number of hits, relative abundance, and coverage.
 * **5.b**: Reduced relative abundance calculation
    * **(5.b.1)** Group recovered MAGs of all samples within the provided metadata table using **gOTUpick**.
    * **(5.b.2)** Select the highest quality MAG within the **gOTUpick** groups as the group's representative MAG.
    * **(5.b.3)** Calculate representative MAGs size, the number of quality-controlled reads, and the average read length of all samples to be mapped.
    * **(5.b.4)** Index selected MAGs
    * **(5.b.5)** Map reads from the sample on the representatives MAGs from the **gOTUpick** group.
    * **(5.b.6)** Calculate representative MAGs absolute number of hits, relative abundance, and coverage tables.
* **5.c**: Genes relative abundance calculation. Currently working only on prokaryotic genes
    * **(5.c.1)** Select MAGs from all samples within the provided metadata table
    * **(5.c.2)** Index MAGs.
    * **(5.c.3)** Map reads from samples on their respective recovered MAGs.
    * **(5.c.4)** Convert the .gff file from **PROKKA** gene annotation on each MAG to .gtf.
    * **(5.c.5)** Calculate average read length of all mapped samples.
    * **(5.c.6)** Calculate gene lenght from **PROKKA** .gtf file.
    * **(5.c.7)** Calculate genes absolute number of hits, relative abudance, coverage, and TPM tables for each MAG.


**UNDER DEVELOPMENT**



# System requirements

**Warning:** complex dependencies ahead!

Hopefully, you won't have to worry too much about it. We designed an installation script that should take care of every dependency for you.
The software used during the pipeline are detailed described here: [Dependencies description](https://github.com/JotaKas/MuDoGeR/blob/master/dependencies_description.md).

Keep in mind that the MuDoGeR pipeline requires some computer power. The complete software installation requires approximately 170 GB, but maker2 uses 99 GB of that space since it requires the database to be installed in a specific manner. See [Module 4 setup](https://github.com/JotaKas/MuDoGeR/blob/master/installation/genemark_maker2_installation.md). The complete database requirements, considering all tools, is around 439.9 GB.

MuDoGeR is designed to support only Linux x64 systems. As for the resource requirements, the MuDoGeR framework uses software that requires a large amount of RAM (e.g **GDTB-Tk**, **MetaWRAP** ). Specific resource requirements vary depending on your data and its sequencing depth. We recommend the user provide at least 180 GB of RAM.
For the assembly process, **MuDoGeR** attempt to calculate the amount of memory necessary for **metaSPades** (on step **1.b**). The user should be aware that samples with higher expected diversity require a higher amount of memory.

# Installation

**1 - Install miniconda**

The MuDoGeR pipeline requires several tools that have multiple dependencies. More often than not, those tools require conflicting dependencies. To tackle this problem, MuDoGeR requires miniconda to be previously installed in your system. MuDoGeR will use miniconda to create multiple environments that are automatically orchestrated during the pipeline's functioning. Following you have a possible way of installing miniconda.

```console

See documentation: https://docs.conda.io/en/latest/miniconda.html

$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

$ chmod +x Miniconda3-latest-Linux-x86_64.sh

$ ./Miniconda3-latest-Linux-x86_64.sh

$ export PATH=~/miniconda3/bin:$PATH

```

**2 - Install MuDoGeR**

Once you have miniconda installed and on your PATH you can properly install MuDoGeR.
The MuDoGeR environment can be complex, and the installation script was designed to install and set up all necessary bioinformatics tools.

```console

#Create the mudoger conda environment
$ conda create mudoger_env

#activate environmet
$ conda activate mudoger_env

#clone repository

$ git clone https://github.com/JotaKas/MuDoGeR.git

$ cd MuDoGeR

#Run installation script as follows
$ bash -i installation/install.sh

#Follow the instructions on the screen

```
**3 - Install necessary databases**

Several bioinformatics tools used within MuDoGeR require specific databases to work. We developed a database download and set up tool to make our lives easier. Make sure to run the database setup after MuDoGeR is installed.

Use this script if you want MuDoGeR to take care of everything. 

```console
#Make sure mudoger_env is activated
$ conda activate mudoger_env

#Go to MuDoGeR cloned directory
$ cd MuDoGeR

#Run database setup script
$ bash -i installation/database-setup.sh -o /path/to/save/databases

```

**4 - Additional module 4 (eukaryotes) installation instructions**

Some tools used in module 4 (GENEMARK and MAKER2) require the user to provide information to the developers. Consequently, we could not implement an automatic installation and setup script. However, we created a tutorial to finish the module 4 setup.

The module 4 setup tutorial is found in [Module 4 setup](https://github.com/JotaKas/MuDoGeR/blob/master/installation/genemark_maker2_installation.md)


# Simplified usage of the MuDoGeR

MuDoGeR was designed to work module by module, starting from pre-process (Module 1). Additional modularity will be added in future updates to allow the user to run specific parts of the pipeline.

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
				Version 1.0.0



Mudoger v=1.0.0
Usage: mudoger --module module_name --meta metadata_table.tsv -o output_folder [module_options]

  --meta              		 Metadata table with your samples, as explained in the github documentation
  --module preprocess            Run all steps from module 1 (read_qc, kmer memory prediction and assembly)
          	read_qc
          	mem_pred
          	assembly
  --module prokaryotes              Recovery of Prokaryotic Metagenome-Assembled Genomes
          	initial_binning
          	bac_bin_ref
          	arc_bin_ref
          	bin_dereplication
          	tax_assignment
          	bin_quality
          	gene_annotation
          	seq_metrics
  --module viruses		Recovery of Uncultivated Viral Genomes
          	uvigs_recov
          	uvigs_tax
          	host_pred
          	uvigs_quality
          	uvigs_metrics
  --module eukaryotes		Recovery of Eukaryotic Metagenome-Assembled Genomes
          	euk_recov
          	gene_pred
          	euk_quality
          	gene_annotation
          	euk_scg
  --module abundance_tables	Abundance calculation
          	brat_type       can be --complete , --reduced or --genes (default: --reduced)
          	mapping_type	can be --absolute-values, --coverage, and --relative-abundance

  --help | -h		show this help message
  --version | -v	show mudoger version

```

**MuDoGeR** is an easy-to-use wrapper of several tools organized within modules. The individual modules can be called independently.

The pipeline requires,as input, a metadata table in tsv format containing the samples to be processed and the path to its raw sequence reads. An example metadata file is as follows:
```
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

A simplified use of MuDoGeR can be done as follows:

```console

$ mudoger --module preprocess --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20

$ mudoger --module prokaryotes ---meta /path/to/metadata.tsv -o /path/to/output/folder -t 20

$ mudoger --module viruses --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20

$ mudoger --module eukaryotes --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20

$ mudoger --module abundance_tables --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20 --coverage --relative-abundance

```

MuDoGeR is designed to completely run all multi-domain genome recovery pipeline. In order for MuDoGeR to work automaticaly, from start to finish, we use a specific folder structure. Please, read the [Manual_MuDoGeR](https://github.com/JotaKas/MuDoGeR/blob/master/Manual_MuDoGeR.md) if you would like to manipulate MuDoGeR. 

The two level folder structure is as follows:
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
    ├── all_bins
    ├── all_metrics
    └── gOTUpick_results
    └── abundance_results

```

A more detailed tutorial for the MuDoGeR can be found in [Manual_MuDoGeR](https://github.com/JotaKas/MuDoGeR/blob/master/Manual_MuDoGeR.md).


# Citing

# Acknowledgements

	



 


	

