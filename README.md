# MuDoGeR

 # Multi-domain Genome Recovery (MuDoGeR)
 
 
 ![](https://github.com/EfthymisF/folder-scripts/blob/master/index.png)

Multi-domain Genome Recovery (MuDoGeR)is a tool developed to help users to recover metagenome assembled genomes from hundreds of samples simultaneously
The MuDoGeR pipeline process combines a number of pipelines from different tools. It starts with **Read Quality control** of the libraries and **Assembly** of the *"good quality reads"* resulted from the former step. For the next step of data analysis, the pipeline is devided in 3 different branches: **Metawrap** pipeline is used for analayzing Prokaryotic genomes. **Virsorter**, **Virfinder** and **Vibrant** pipelines are followed and combined for the prediction of viral sequences. **EukRep** pipeline is run for metagenomic analysis of Eukaryotic genomes. 
 

## MuDoGeR Overview


### (1) Pre-Processing
* **(1.1)** Trimming of the metagenomic library and removal of all the host reads by running  metaWRAP-Read_qc
* **(1.2)** Assembly of the clean metagenomic reads produced in **(1.1)** by running metaWRAP-Assembly module

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

###  (4) Metagenomic recovery of Eukaryotic genomes
* **(4.1)** Classification  of  Eukaryotic genomes with EukRep.
* **(4.2)** Use of CONCOCT for binning and filtering the Eukaryotic sequences by size.
* **(4.3)** Then, genes are predicted by using the trained GeneMark-ES model 
* **(4.4)** MAKER2. 
* **(4.5)** BUSCO is applied for detection of single copy orthologous genes (SCGs) and will estimate the completeness and contamination of Eukaryotic reads



# Parts of the MuDoGeR

* **Installation** 
* **Read Quality control:** Trimming of the reads and removal of possible host reads
* **Assembly:** Assembly of "good quality" sequences
* **MetaWRAP pipeline:**  Metagenomic analysis of prokaryotic genomes
* **Pipelines for viral genomes(VirFinder, VirSorter, Vibrant):** Metagenomic analysis of viral sequences 
* **EukRep pipeline:** Metagenomic analysis of Eukaryotic genomes


# System requirements


# Installation

### MetaWRAP

The instructions for the metaWrap installation can be found at ![metaWRAP_Manual](https://github.com/EfthymisF/new/blob/master/Tutorial.md)

### VirFinder

The instructions for Virfinder installation can be found at ![VirFinder_Manual](https://github.com/jessieren/VirFinder)

### VirSorter

The instructions for VirSorter installation can be found at ![VirSorter_Manual](https://github.com/simroux/VirSorter) 

### VIBRANT 

The instructions for VIBRANT installation can be found at ![VIBRANT_Manual](https://github.com/AnantharamanLab/VIBRANT) 

### EukRep
The instructions for EukRep installation can be found at ![EukRep_Manual](https://github.com/patrickwest/EukRep_Pipeline ) 



# Using MuDoGeR

A more detailed tutorial for the MuDoGeR can be found in ![Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md) file.Once all the dependencies are in place, running metaWRAP is relatively simple. The main metaWRAP script wraps around all of its individual modules, which you can call independently. The help message for MuDoGeR is the following:

```mudoger -h
	Usage: mudoger [module] --help
	Options:
	
	Full_Workflow	Whole workflow excecution
	read_qc		Raw read QC module
	assembly	Assembly module
	Prokaryotic	Recovery of Prokaryotic genomes using metaWRAP
	Viral		Recovery of Viral genomes using VirFinder, VirSorter, VIBRANT
	Eukaryotic 	Recovery of Eukaryotic genomes using EukRep
	Sep_Modules	Watch each module separately
	
	--help | -h		show this help message
	--version | -v	show MuDoGeR version
	--show-config	show where the mudoger configuration files are stored

	
	
 ```


# Step 0:Dowloading the libraries

# Step 1:Pre-Processing 
## Step 1.1: Read Quality control
## Step1.2: Assembly of the good quality control reads

# Step 2:Metagenomic recovery of Prokaryotic genomes 

# Step 3:Pipelines for viral genomes recovery(VirFinder, VirSorter, Vibrant)

The help message of the viral recovery script is the following: 


```viral -h

	Usage: Usage: Viral assembly and dereplication [options] -1 P19011_A1_assembly.fasta -o output_dir --help
	Options:
	
	-f STR          Assembly.fasta
	-o STR          output directory
	
	-virfinder     Recovery of viral data with VirFinder
	-virsorter     Recovery of viral data with VirSorter
	-vibrant       Recovery of viral data with VIBRANT
	-dereplication Removal of replicate sequences
	
	--help | -h		show this help message
	--version | -v	show virfinder,virsorter,vibrant version
	--show-config	show where the mudoger configuration files are stored

	
	
 ```

* The recovery of viral genomes is achieved by using the tolls VirFinder, VirSorter and VIBRANT. Although each of the tools can be used independently, to achieve the maximum recovery of viral genomes it is better to use all the three tools together. After the genome recovery the, the viral genomes recovered from each method are filtered:
* In the VirFinder results sequences with q-value > 0.01 and lentgh =< 1000 bp are removed
* In the VirSorter results only the sequences of the categories 1 and 2 are chosen 
* In the VIBRANT results only the sequences of the `Assembly.phages_combined.fna` are chosen

* The headers of the filtered sequences are combined to one file and the repeated sequences are removed. With the help of the headers, those sequences are extracted from the initial Assembly.fasta file and are located to a new fasta file. Finally the replicates are removed with a dereplication process. The loading parameters of this part of the step are maximum coverage `-c=70` and maximal identity `-i=95`. The combination of identity and coverage parameters is able to vary, depending on the purposes and the goals of the user. 

* Run the pipelines for viral genomes recovery module

* ``` ~/viral_working_directory/viral_recovery_script -f ~/path/to/assembly/file -o /path/to/output/folder -o /path/to/output/folder/file/virsorter.tsv``` 





 
# Step 4:EukRep pipeline



 


	

