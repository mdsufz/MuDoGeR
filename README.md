# MuDoGeR

 # Multi-domain Genome Recovery (MuDoGeR)
 
 
 ![](https://github.com/mdsufz/MuDoGeR/blob/master/Mudoger_scheme_update.png)

The **Multi-domain Genome Recovery (MuDoGeR)** pipeline is a tool developed to help users to recover Metagenome-Assembled Genomes (MAGs (Parks, 2018)) from dozens to hundreds of samples simultaneously. The **MuDoGeR** pipeline combines pipelines from different tools. The **MuDoGeR** pipeline starts with **Pre-Processing**, which is divided in 2 steps: The first one, **Read Quality Control**, feeds the **Assembly** and the assemblies will be used for the other pipelines. After pre-processing of the data, **MuDoGeR** is divided in 3 different branches: **Recovery of Prokaryotic Metagenome-Assembled Genomes** (using **MetaWrap**), **Recovery of Viral Metagenome-Assembled Genomes** (using **VirSorter**, **VirFinder** and **VIBRANT**) and **Recovery of Eukaryotic Metagenome-Assembled Genomes** (using **EukRep**). Furthermore, a strategy was developed for mapping the relative abundance of the selected maps in each library. Also, a step was added for bin_refinement of the selected Metagenome-Assembled Genomes from Prokaryotes, using **U-bin**.
 

## MuDoGeR Overview

* Instructions for using the **MuDoGeR** can be found in the following hyperlink: ![Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md).
* Detailed description of the **MuDoGeR** steps can be found in the following hyperlink: ![Module description](https://github.com/mdsufz/MuDoGeR/blob/master/module_description.md).
* Detailed instructions for the installation of the **MuDoGeR** tools can be found in the following hyperlink: ![Installation](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#installation).
* Information about the system requirements of the **MuDoGeR** can be found in the following hyperlink: ![System requirements](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#system-requirements).
* The simplified usage of the **MuDoGeR** can be found in the following hyperlink: ![Simplified usage of the MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#using-mudoger).

### (1) Pre-Processing

![](https://github.com/mdsufz/MuDoGeR/blob/master/Preprocess_module_update.png)

**Figure 2.** Module 1 of the MuDoGeR pipeline.

 The steps of Module 1 are shown in Figure 2 and can be excecuted with the scripts find in the following hyperlink: ![Pre-Processing module](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-1-pre-processing-module)
 
 The steps of Module 1 can be summarized as following:

* **(1.a)** Trimming of the metagenomic library and removal of all the host reads by running  **metaWRAP-Read_qc**.
* **(1.b)** Calculation of the amount of resources for the good quality reads of the libraries.
	* **(1.b.1)** The k-mer of the good quality reads produced in **(1.a)** is calculated. The k-mer sizes that will be investigated are 33 and 55.  
	* **(1.b.2)** The calculated k-mer is added to an equation that is used to estimate the amount of memory that **MetaSpades** utilizes to assemble the good quality reads.
* **(1.c)** Assembly of the good quality reads with **MetaSpades**.

### (2) Recovery of Prokaryotic Metagenome-Assembled Genomes

The different steps of the module 2 are shown in Figure 3 and excecuted with the scripts find in the following hyperlink: ![Pipeline for recovery of Prokaryotic Metagenome-Assembled Genomes](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-2-pipelines-for-prokaryotic-genome-recovery).

 The steps of Module 2 can be summarized as following:

* **(2.a)** Binning with **MaxBin2**, **metaBAT2**, and **CONCOCT**. 
* **(2.b)** Dereplication of bins for prior bin_refinement. The completeness/contamination parameters have been set to 50%/10% for Bacteria and 40%/30% for Archaea. 
* **(2.c)** Taxonomic classification of the genomic bins produced in **(2.b)** using **GTDB-Tk**.
* **(2.d)** Generation of quality matrix of genomic bins produced in **(2.b)** using **CheckM**. 
* **(2.e)** Filtering of genomic bins produced in **(2.b)** by quality (Parks (2018)). The minimum quality for the filtering is set by default at 50. In this step, the user can also change the required quality (optional step).      
* **(2.f)** Annotation of genomic bins produced in **(2.b)** with **PROKKA**.
* **(2.g)** Selection of Prokaryotic Metagenome-Assembled Genomes.
* **(2.h)** Refinement of Prokaryotic Metagenome-Assembled Genomes using **U-bin** (optional step).

###  (3) Recovery of Viral Metagenome-Assembled Genomes

![](https://github.com/mdsufz/MuDoGeR/blob/master/VIral_pipeline_.png)
**Figure 4.** Module 3 of the MuDoGeR pipeline.

The steps of the module 3  are shown in Figure 4 and excecuted in the scripts find in the following hyperlink: ![Pipelines for viral genomes recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-3-pipelines-for-viral-genomes-recovery)

 The steps of Module 3 can be summarized as following:

* **(3.a)** Recovery of viral genome metagenomes using **Virfinder**, **VirSorter** and **VIBRANT**.
* **(3.b)** Filtering of the recovered genomes.
* **(3.c)** Combination of the headers of all the Viral Metagenome-Assembled Genomes in a single file. Then, removal of the headers of the repeated sequences and sorting the remaining headers by sequence length.
* **(3.e)** Assembly of the headers in the file produced in **(3.c)** with the respective sequences from the assembly file produced by **MetaSpades** in **(1.3)**, for the generation of a fasta file.
* **(3.f)** Removal of replicates from the assemblies in the extracted file using de-replication function with minimum coverage of 70% and minimum identity of 95%.
* **(3.g)** Checking the quality of the dereplicated contigs with **CheckV**.
* **(3.h)** Taxonomic classification of the dereplicated contigs with **vContact2**. 
* **(3.i)** Host identification of the dereplicated contigs using **WIsH**. 
* **(3.j)** Selection of Viral Metagenome-Assembled Genomes.

###  (4) Recovery of Eukaryotic Metagenome-Assembled Genomes

![](https://github.com/mdsufz/MuDoGeR/blob/master/Eykaryotic_module.png)
**Figure 5.**  Module 4 of the **MuDoGeR** pipeline.

The steps of the module 4  are shown in Figure 5 and can be excecuted with the scripts found in ![Pipelines for eukaryotic genomes recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-4-pipelines-for-eukaryotic-genomes-recovery)

 The steps of Module 4 can be summarized as following:

* **(4.a)** Classification of Eukaryotic assemblies and removal of prokaryotic assemblies with **EukRep**.
* **(4.b)** Use of **CONCOCT** for binning the Eukaryotic assemblies.  
* **(4.c)** Filtering the Eukaryotic bins, produced by **CONCOCT**, by size. Bins with size <= 2.5 Mb are removed.
* **(4.d)** In the filtered bins produced in **(4.c)**, genes are predicted by the trained **GeneMark-EV** model.
* **(4.e)** **EukCC** utilization for estimating the contamination of Eukaryotic filtered bins produced in **(4.c)**.
* **(4.f)** **MAKER2** annotates the predicted genes produced by **GeneMark-EV**. 
* **(4.g)** **BUSCO** is applied to the annotated genes from **MAKER2**, for detection of single copy orthologous genes (SCGs) and estimation of the completeness of Eukaryotic contigs.
* **(4.h)** Selection of Eukaryotic Metagenome-Assembled Genomes.

# System requirements
(not done yet)

# Installation

## Prokaryotic module
### MetaWRAP
The instructions for the metaWrap installation can be found at ![metaWRAP_Manual](https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md).
### uBin
The instructions for the uBin installation can be found at ![uBin_Manual](https://github.com/ProbstLab/uBin-helperscripts.).

## Viral module

### VirFinder
The instructions for Virfinder installation can be found at ![VirFinder_Manual](https://github.com/jessieren/VirFinder).
### VirSorter
The instructions for VirSorter installation can be found at ![VirSorter_Manual](https://github.com/simroux/VirSorter). 
### VIBRANT 
The instructions for VIBRANT installation can be found at ![VIBRANT_Manual](https://github.com/AnantharamanLab/VIBRANT). 
### vContact
The instructions for vContact installation can be found at ![vContact_Manual](https://github.com/geeklhem/vcontact/blob/master/doc/source/tutorial.rst).

## Eukaryotic module
### EukRep
The instructions for EukRep installation can be found at ![EukRep_Manual](https://github.com/patrickwest/EukRep_Pipeline ). 
### EukCC
The instructions for EukCC installation tool can be found at ![EukCC_Manual](https://github.com/Finn-Lab/EukCC).



# DATABASES 
(Not done yet)


|   Database	|     Used 	|
| ------------- | ------------- |
|   Checkm_DB	| Prokaryotic module |
| KRAKEN standard database | Prokaryotic module  |
| First Header  | Prokaryotic module  |
| KRAKEN2 standard database  |Prokaryotic module |
| NCBI_nt  | Prokaryotic module |
| NCBI_tax  | Prokarytic module |
| Indexed hg38 | Read_qc module |



# Simplified usage of the MuDoGeR

A more detailed tutorial for the MuDoGeR can be found in ![Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md) file. In this file, instructions and examples can be found.

After all dependencies are in position then the MuDoGeR is ready to run. The run of MuDoGeR is easy and there is a main script that can wrap around all the MuDoGeR's individual module. The individual modules can be called independently.

```
mudoger -h
	Usage: mudoger [module] --help
	Options:
	
	read_qc		Raw read QC module
	assembly	Assembly module
	resources	Calculation of memory required by the MetaSpades, for assembling the good quality reads 
	prokaryotic	Recovery of Prokaryotic genomes using metaWRAP module
	viral		Recovery of Viral genomes using VirFinder, VirSorter, VIBRANT module
	eukaryotic 	Recovery of Eukaryotic genomes using EukRep module
		
	--help | -h		show this help message
	--version | -v	show MuDoGeR version
	--show-config	show where the mudoger configuration files are stored

 ```

Each of the modules is run seperately. As an example to run the viral module:

```
mudoger viral module -h

Usage: mudoger viral module [options] -o output_dir -f assembly.fasta 
Options:
	
	-o STR          output directory
	-f STR          assembly fasta file
	-c INT		minimum coverage (default=70)
	-i INT		minimum identity (default=95)
	
	-virfinder	Recovery of viral data with VirFinder 
	-virsorter	Recovery of viral data with VirSorter 	
	-vibrant	Recovery of viral data with VIBRANT
	-dereplication	Removal of replicate sequences
	-checkv		Quality control of dereplicated contigs with CheckV		
	-vcontact2	Taxonomy of dereplicated contigs with vContact2
	--wish		Host identification of dereplicated contigs with WIsH
```

# Citing

# Acknowledgements

	



 


	

