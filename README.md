# MuDoGeR

 # Multi-domain Genome Recovery (MuDoGeR)
 
 
 ![](https://github.com/mdsufz/MuDoGeR/blob/master/mudoger.png)

Multi-domain Genome Recovery (MuDoGeR)is a tool developed to help users to recover metagenome assembled genomes from hundreds of samples simultaneously
The MuDoGeR pipeline process combines a number of pipelines from different tools. It starts with **Read Quality control** of the libraries and **Assembly** of the *"good quality reads"* resulted from the former step. For the next step of data analysis, the pipeline is devided in 3 different branches: **Metawrap** pipeline is used for analayzing Prokaryotic genomes. **Virsorter**, **Virfinder** and **Vibrant** pipelines are followed and combined for the prediction of viral sequences. **EukRep** pipeline is run for metagenomic analysis of Eukaryotic genomes. 
 

## MuDoGeR Overview

For the usage of MuDoGeR the user can follow the instructions in the ![Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md) file
For detailed description of the MuDoGeR steps, the user can study the ![Module description](https://github.com/mdsufz/MuDoGeR/blob/master/module_description.md) file

### (1) Pre-Processing
* **(1.1)** Trimming of the metagenomic library and removal of all the host reads by running  metaWRAP-Read_qc
* **(1.2)** Kmer calculation of the clean metagenomic reads produced in **(1.1)** module **(1.2.1)** and assembly of the good quality reads by running metaWRAP-Assembly module **(1.2.2)**

###  (2) Metagenomic recovery of Prokaryotic genomes 

The steps of the module are shown in Figure 3 and can be excecuted with the scripts found in ![Pipelines for prokaryotic genome recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-2-pipelines-for-prokaryotic-genome-recovery)


* **(2.1)** Bin extraction with MaxBin2, metaBAT2, and CONCOCT 
* **(2.2)** Unification of multiple binning prediction into an advanced bin set (Bin refinement) 
* **(2.3)** Quality control using CheckM 
* **(2.4)**
* **(2.5)** Classification of genomic bins with GTDB
* **(2.6)** Reassemble the final annotation with PROKKA
* **(2.7)** Selection of MAGs
* **(2.8)** U-bin tool for refinement of prokaryotes




###  (3) Metagenomic recovery of Viral genomes
The steps of the module 3  are shown in Figure 4 and can be excecuted with the scripts found in ![Pipelines for viral genomes recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-3-pipelines-for-viral-genomes-recovery)

![](https://github.com/mdsufz/MuDoGeR/blob/master/Viral%20Module.png)

* **(3.1)** Recovery of viral metagenomes using Virfinder, VirSorter and VIBRANT for the prediction of viral sequences . 
* **(3.2)** Filtering of the recovered genomes
* **(3.3)** Combination of the headers of the filtered data to a single file, removed repeated sequences and sort by length
* **(3.4)** Extraction of the viral sequences from the assembly file using the headers
* **(3.5)** Removal of replicates with de-replication function
* **(3.6)** Checking the quality of the dereplicated contigs with CheckV
* **(3.7)** Taxonomy of of the dereplicated contigs with vContact
* **(3.8)** Host identification of the dereplicated contigs using WIsH 


###  (4) Metagenomic recovery of Eukaryotic genomes
* **(4.1)** Classification  of  Eukaryotic genomes with EukRep.
* **(4.2)** Use of CONCOCT for binning and filtering the Eukaryotic sequences by size. The filterin
* **(4.3)** Then, genes are predicted by using the trained GeneMark-ES model 
* **(4.4)** MAKER2 for annotation
* **(4.5)** BUSCO is applied for detection of single copy orthologous genes (SCGs) and will estimate the completeness of Eukaryotic contigs
* **(4.6)** EukCC for estimating the contamination of Eukaryotic contigs


# Parts of the MuDoGeR

* **Installation** 
* **Read Quality control:** Trimming of the reads and removal of possible host reads
* **Assembly:** Assembly of "good quality" sequences
* **MetaWRAP pipeline:**  Metagenomic analysis of prokaryotic genomes
* **Pipelines for viral genomes(VirFinder, VirSorter, Vibrant):** Metagenomic analysis of viral sequences 
* **EukRep pipeline:** Metagenomic analysis of Eukaryotic genomes



# System requirements


# Installation

## Prokaryotic module
### MetaWRAP

The instructions for the metaWrap installation can be found at ![metaWRAP_Manual](https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md)

## Viral module

### VirFinder
The instructions for Virfinder installation can be found at ![VirFinder_Manual](https://github.com/jessieren/VirFinder)
### VirSorter
The instructions for VirSorter installation can be found at ![VirSorter_Manual](https://github.com/simroux/VirSorter) 
### VIBRANT 
The instructions for VIBRANT installation can be found at ![VIBRANT_Manual](https://github.com/AnantharamanLab/VIBRANT) 

## Eukaryotic module
### EukRep
The instructions for EukRep installation can be found at ![EukRep_Manual](https://github.com/patrickwest/EukRep_Pipeline ) 



# DATABASES


|   Database	|     Used 	|
| ------------- | ------------- |
|   Checkm_DB	| Prokaryotic module |
| KRAKEN standard database | Prokaryotic module  |
| First Header  | Prokaryotic module  |
| KRAKEN2 standard database  |Prokaryotic module |
| NCBI_nt  | Prokaryotic module |
| NCBI_tax  | Prokarytic module |
| Indexed hg38 | Read_qc module |



# DETAILED PIPELINE WALKTHROUGH

# Using MuDoGeR

A more detailed tutorial for the MuDoGeR can be found in ![Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md) file. In this file, instructions and examples can be found.

After all dependencies are in position then the MuDoGeR is ready to run. The run of MuDoGeR is easy and there is a main script that can wrap around all the MuDoGeR's individual module. The individual modules can be called independently.

```
mudoger -h
	Usage: mudoger [module] --help
	Options:
	
	read_qc		Raw read QC module
	assembly	Assembly module
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
	
	--virfinder	Recovery of viral data with VirFinder 
	--virsorter	Recovery of viral data with VirSorter 	
	--vibrant	Recovery of viral data with VIBRANT
	--dereplication	Removal of replicate sequences
	--checkv	Quality control of dereplicated contigs with CheckV		
	--vcontact2	Taxonomy of dereplicated contigs with vContact2
	--wish		Host identification of dereplicated contigs with WIsH
```

# Citing

# Acknowledgements

	



 


	

