# MuDoGeR

 # Multi-domain Genome Recovery (MuDoGeR)
 
 
 ![](https://github.com/mdsufz/MuDoGeR/blob/master/Mudoger_scheme_update.png)

Multi-domain Genome Recovery (MuDoGeR)is a tool developed to help users to recover metagenome assembled genomes from hundreds of samples simultaneously
The MuDoGeR pipeline process combines a number of pipelines from different tools. It starts with **Read Quality control** of the libraries and **Assembly** of the *"good quality reads"* resulted from the former step. For the next step of data analysis, the pipeline is devided in 3 different branches: **Metawrap** pipeline is used for analayzing Prokaryotic genomes. **Virsorter**, **Virfinder** and **Vibrant** pipelines are followed and combined for the prediction of viral sequences. **EukRep** pipeline is run for metagenomic analysis of Eukaryotic genomes. Also, a strategy was developed for mapping the relative abundance of the selected maps in each library.
 

## MuDoGeR Overview

* For the usage of MuDoGeR the user can follow the instructions in the ![Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md) file
* For detailed description of the MuDoGeR steps, the user can study the ![Module description](https://github.com/mdsufz/MuDoGeR/blob/master/module_description.md) file

### (1) Pre-Processing


![](https://github.com/mdsufz/MuDoGeR/blob/master/Pre_processing.png)

Figure 2. **Pre-Processing** 

The steps of the module 21 are shown in Figure 2 and can be excecuted with the scripts found in: ![Pre-Processing module](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-1-pre-processing-module)

* **(1.a)** Trimming of the metagenomic library and removal of all the host reads by running  metaWRAP-Read_qc
* **(1.b)** Calculation of the amount of resources for the good quality reads of the libraries.
* **(1.b.1)** The k-mer of the good quality reads produced in **(1.1)** is first calculated.The k-mer size which will investigated is 33 and 55.
* **(1.b.2)** The calculated k-mer is added to an equation which will estimate the amount of memory that metaSpades will utilize to assemble the good quality reads. 
* **(1.c)** Assembly of the good quality reads with metaSpades

###  (2) Metagenomic recovery of Prokaryotic genomes 

The steps of the module 2 are shown in Figure 3 and can be excecuted with the scripts found in ![Pipelines for prokaryotic genome recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-2-pipelines-for-prokaryotic-genome-recovery)


* **(2.a)** Bin extraction with MaxBin2, metaBAT2, and CONCOCT 
* **(2.b)** Unification of multiple binning prediction into an advanced bin set (Bin refinement). The completeness/contamination parameters have been setted to 50/10 for Bacteria and 40/30 for archaea 
* **(2.c)** Classification of genomic bins produced in **(2.b)** (or **(2.d)** if the user chooses to) step with GTDB.
* **(2.d)** Quality control of the bins produced in **(2.b)** using CheckM 
* **(2.e)** Filtering of the bins produced in **(2.b)** with higher completeness than 50. The user here can also change the size of the required completeness (optional step)      
* **(2.f)** Reassemble the final annotation of the bins produced in  **(2.b)** with PROKKA
* **(2.g)** Selection of prokaryotic MAGs
* **(2.h)** U-bin tool for manual curation of of genomes form the chosen MAGs 




###  (3) Metagenomic recovery of Viral genomes



Figure 4. **Viral pipeline**

The steps of the module 3  are shown in Figure 4 and can be excecuted with the scripts found in ![Pipelines for viral genomes recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-3-pipelines-for-viral-genomes-recovery)

* **(3.a)** Recovery of viral metagenomes using Virfinder, VirSorter and VIBRANT for the prediction of viral sequences . 
* **(3.b)** Filtering of the recovered genomes
* **(3.c)** Combination of the headers of the filtered data to a single file, removed repeated sequences and sort by length
* **(3.e)** Extraction of the viral sequences from the assembly file using the headers
* **(3.f)** Removal of replicates from the assemblies in the extracted file using de-replication function
* **(3.g)** Checking the quality of the dereplicated contigs with CheckV
* **(3.h)** Taxonomy of the dereplicated contigs with vContact 
* **(3.i)** Host identification of the dereplicated contigs using WIsH 
* **(3.j)** Selection of viral MAGs


###  (4) Metagenomic recovery of Eukaryotic genomes

The steps of the module 4  are shown in Figure 5 and can be excecuted with the scripts found in ![Pipelines for eukaryotic genomes recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-4-pipelines-for-eukaryotic-genomes-recovery)

* **(4.a)** Classification of Eukaryotic assemblies and removal of prokaryotic assemblies with EukRep
* **(4.b)** Use of CONCOCT for binning the Eukaryotic assemblies  
* **(4.c)** Filtering the Eukaryotic bins produced by CONCOCT by size. Bins with size <= 2.5 Mb are removed
* **(4.d)** In the filtered bins produced in **(4.c)**, genes are predicted by the trained GeneMark-ES model   
* **(4.e)** MAKER2 annotates the predicted genes produced by GeneMark-ES 
* **(4.f)** BUSCO is applied to the annotated genes from for MAKER2, for detection of single copy orthologous genes (SCGs) and estimation of the completeness of Eukaryotic contigs
* **(4.g)** EukCC utilization for estimating the contamination of Eukaryotic filtered bins produced in **(4.c)** 
* **(4.h)** Selection of eukaryotic MAGs


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
### uBin
The instructions for the uBin installation can be found at ![uBin_Manual](https://github.com/ProbstLab/uBin-helperscripts.)

## Viral module

### VirFinder
The instructions for Virfinder installation can be found at ![VirFinder_Manual](https://github.com/jessieren/VirFinder)
### VirSorter
The instructions for VirSorter installation can be found at ![VirSorter_Manual](https://github.com/simroux/VirSorter) 
### VIBRANT 
The instructions for VIBRANT installation can be found at ![VIBRANT_Manual](https://github.com/AnantharamanLab/VIBRANT) 
### vContact
The instructions for vContact installation can be found at ![vContact_Manual](https://github.com/geeklhem/vcontact/blob/master/doc/source/tutorial.rst)

## Eukaryotic module
### EukRep
The instructions for EukRep installation can be found at ![EukRep_Manual](https://github.com/patrickwest/EukRep_Pipeline ) 
### EukCC
The instructions for EukCC installation tool can be found at ![EukCC_Manual](https://github.com/Finn-Lab/EukCC)



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

	



 


	

