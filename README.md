## thym


 # Multidomain Genome Recovery (MUDOGER)
 
 
 ![](https://github.com/EfthymisF/folder-scripts/blob/master/index.png)
 
 MUDOGER is a tool developed for metagenomic analysis and classification of data contained in mixed libraries from a wide range of environments. MUDOGER is user-friendly  and is the first tool that combines a number of pipelines from diffrent tools for the simultaneous extraction of  Prokaryotic, Eykaryotic and Viral data. The MUDOGER pipeline process starts with **Read Quality control** of the libraries and **Assembly** of the *"good quality reads"* resulted from the former step. For the next step of  data analysis, the pipeline is devided in 3 different branches: **1. Metawrap** pipeline is used for analayzing Prokaryotic genomes. **2.** **Virsorter**, **Virfinder** and **Vibrant** pipelines are followed and combined for the prediction of viral sequences. **3.** **Eykarep** is run for metagenomic analysis of Eykaryotic genomes. 
 

## MUDOGER Workflow

##  (1) Trimming of the Metagenomic library and removal of all the host reads by running  metaWRAP-Read_qc

##  (2) Assembly of the clean metagenomic reads produced in **(1)** by running metaWRAP-Assembly module

##  (3) Metagenomic analysis of Prokaryotic genomes 
```(3.1) Bin extraction with MaxBin2, metaBAT2, and CONCOCT (3.2) Unification of multiple binning prediction into an advanced bin set (Bin refinement) (3.3) Quality control using cCheckM (3.4) Classification of genomic bines (3.5) Reassemblethe final annotation with PROKKA```

##  (4) Metagenomic analysis of Viral genomes
```(4.1) Recovery of viral metagenomes using Virfinder, Virsorter and Vibrant for the prediction of viral sequences, combination of them to a single file and removal of replicates (4.2) Prediction of each contig's protein by Prodigal (4.3) Identification of Proteins by Blastp``` 

##  (5) Metagenomic analysis of Eykaryotic genomes
```(5.1) Classification  of  Eykaryotic genomes with EykRep (5.2) Use of CONCOCT for binning and filtering of the Eykaryotic sequences by size.(5.3) Then, genes are predicted by  using the trained GeneMark-ES model and MAKER2. (5.4) BUSCO is applied for detection of single copy orthologous genes (SCGs) and will estimate the completeness and contamination of Eykaryotic reads```

             


# OVERVIEW OF THE MUDOGER

* **Installation** 
* **ReadQC Quality control:** Trimming of the reads and removal of possible host reads
* **Assembly:** Assembly of "good quality" sequences
* **METAWRAP pipeline:**  Metagenomic analysis of prokaryotic genomes
* **Pipelines for viral genomes(VirFinder, VirSorter, Vibrant):** Metagenomic analysis of viral sequences 
* **EykRep pipeline:** Metagenomic analysis of Eykaryotic genomes


# System requirements


# Installation

## METAWRAP
```  
* conda create -y -n metawrap-env python=2.7
* conda activate metawrap-env

* # Note: ordering is important 
* conda config --add channels defaults
* conda config --add channels conda-forge
* conda config --add channels bioconda
* conda config --add channels ursky
```




# Using MUDOGER
A tutorial of MUDOGER usage can be found in  ![Manual](https://github.com/EfthymisF/new/blob/master/Tutorial.md)

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
	dereplicate	Removal of redundant viral particles module
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

# Step 3:METAWRAP pipeline

# Step 4:Pipelines for viral genomes(VirFinder, VirSorter, Vibrant)

# Step 5:EykRep pipeline 


 


	

