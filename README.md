# MuDoGeR

 # Multi-domain Genome Recovery (MuDoGeR)
 
 
 ![](https://github.com/EfthymisF/folder-scripts/blob/master/index.png)

Multi-domain Genome Recovery (MuDoGeR)is a tool developed to help users to recover metagenome assembled genomes from hundreds of samples simultaneously
The MuDoGeR pipeline process combines a number of pipelines from different tools. It starts with **Read Quality control** of the libraries and **Assembly** of the *"good quality reads"* resulted from the former step. For the next step of data analysis, the pipeline is devided in 3 different branches: **Metawrap** pipeline is used for analayzing Prokaryotic genomes. **Virsorter**, **Virfinder** and **Vibrant** pipelines are followed and combined for the prediction of viral sequences. **EykRep** pipeline is run for metagenomic analysis of Eykaryotic genomes. 
 

## MuDoGeR Overview


### (1) Pre-Processing
* **(1.1)** Trimming of the Metagenomic library and removal of all the host reads by running  metaWRAP-Read_qc
* **(1.2)** Assembly of the clean metagenomic reads produced in **(1)** by running metaWRAP-Assembly module

###  (2) Metagenomic recovery of Prokaryotic genomes 
* **(2.1)** Bin extraction with MaxBin2, metaBAT2, and CONCOCT 
* **(2.2)** Unification of multiple binning prediction into an advanced bin set (Bin refinement) 
* **(2.3)** Quality control using cCheckM 
* **(2.4)** Classification of genomic bins 
* **(2.5)** Reassemble the final annotation with PROKKA

###  (3) Metagenomic recovary of Viral genomes
* **(3.1)** Recovery of viral metagenomes using Virfinder, VirSorter and VIBRANT for the prediction of viral sequences, combination of them to a single file and removal of replicates. 
* **(3.2)** Prediction of each contig's protein using Prodigal. 
* **(3.3)** Identification using Proteins by Blastp. 

###  (4) Metagenomic recovery of Eykaryotic genomes
* **(4.1)** Classification  of  Eykaryotic genomes with EykRep.
* **(4.2)** Use of CONCOCT for binning and filtering the Eykaryotic sequences by size.
* **(4.3)** Then, genes are predicted by using the trained GeneMark-ES model 
* **(4.4)** MAKER2. 
* **(4.5)** BUSCO is applied for detection of single copy orthologous genes (SCGs) and will estimate the completeness and contamination of Eykaryotic reads



# OVERVIEW OF THE MUDOGER

* **Installation** 
* **Read Quality control:** Trimming of the reads and removal of possible host reads
* **Assembly:** Assembly of "good quality" sequences
* **MetaWRAP pipeline:**  Metagenomic analysis of prokaryotic genomes
* **Pipelines for viral genomes(VirFinder, VirSorter, Vibrant):** Metagenomic analysis of viral sequences 
* **EykRep pipeline:** Metagenomic analysis of Eykaryotic genomes


# System requirements


# Installation

### MetaWRAP

The instructions for the metaWrap installation can be found at ![MetaWrap_Manual](https://github.com/EfthymisF/new/blob/master/Tutorial.md)

### VirFinder

The instructions for Virfinder installation can be found at ![VirFinder_Manual](https://github.com/jessieren/VirFinder)

### VirSorter

The instructions for VirSorter installation can be found at ![VirSorter_Manual](https://github.com/simroux/VirSorter) 

### VIBRANT 

The instructions for VIBRANT installation can be found at ![VIBRANT_Manual](https://github.com/AnantharamanLab/VIBRANT) 

### EykRep
The instructions for EykRep installation can be found at ![EykRep_Manual](https://github.com/patrickwest/EukRep_Pipeline ) 



# Using MuDoGeR

A tutorial of MUDOGER usage can be found in ![Manual](https://github.com/EfthymisF/new/blob/master/Tutorial.md)

Once all the dependencies are in place, running metaWRAP is relatively simple. The main metaWRAP script wraps around all of its individual modules, which you can call independently.

```MUDOGER -h
	Usage: MUDOGER [module] --help
	Options:

	read_qc		Raw read QC module
	assembly	Assembly module
 	binning		Binning of Prokaryotic data module 
	bin_refinement	Refinement of bins from binning module
	checkm		Completeness and contamination of Prokaryotic bins module
	classify 	Classification of  Prokaryotic bins module
	prokka		Annotation of bins module
 	virfinder   	Prediction of viral genomes with VirFinder module 
	virsorter	Prediction of viral genomes with VirSorter module 
	vibrant		Prediction of viral genomes with VIBRANT module
 	viral_filt	Selecion of virfinder,virsorter,vibrant good contings 
	viral_gc_comb	Combination of good viral contigs to one database module	
	dereplication	Removal of redundant viral particles module
	prodigal	Prodigal module	
	blastp		Protein identification module
	eykrep		EykRep classification module
	concoct		CONCOCT module
	genemark	GeneMark-ES module
	maker		Maker2 module
	busco		BUSCO module
	
 ```


# Step 0: Dowloading the libraries

# Step 1: Trimming the reads and removing possible host reads

# Step 2: Assembly of "good quality" sequences

# Step 3:MetaWRAP pipeline

# Step 4:Pipelines for viral genomes(VirFinder, VirSorter, Vibrant)

# Step 5:EykRep pipeline 


 


	

