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
* **2.b**: Taxonomic classification, quality estimation, selection of Prokaryotic MAGs.
     * **(2.b.1)** Taxonomic classification of the prokaryotic bins produced in **(2.a.3)** using **GTDB-Tk**.
     * **(2.b.2)** Generation of quality matrix of the prokaryotic bins produced in **(2.a.3)** using **CheckM**. 
     * **(2.b.3)** Selection of the prokaryotic MAGs using the information retrieved in **(2.b.1)** and **(2.b.2)**.
* **2.c**: Sequence metrics calculation and gene annotation.
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
 * **3.b**: Taxonomic annotation and Quality estimation of Uvigs
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
* **4.c**: Selection of eMAGS. (CHECK HERE IF NECESSARY!!!!!!!!!!!!!!!)
	* **(4.c.1)** Grouping of eMAGs using **BUSCO** and **EukCC** results.
	* **(4.c.2)** Separation of the cluster produced in **(4.c.1)** by species, using ANI (Average Nucleotide Identity) splitter, with deafault ANI_distance 0.95.
	* **(4.c.3)** Selection of representative Metagenome-Assembled Genomes, from the clusters produced in **(4.c.2)**.

### Module 5 Relative abundace
[](https://github.com/JotaKas/MuDoGeR/blob/master/Module%205.PNG)

The steps of the Module 5 are shown in Figure 6. A detailed description of its execution and outputs are found here: [Pipelines for Relative abundace ](https://github.com/JotaKas/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-5-relative-abundance).
* **5.a** Calculation of relative abundance and genome coverage of Prokaryotic Metagenome-Assembled Genomes and construction of relative abundance and genome coverage tables
	* **(5.a.1)** Merging of Paired-End (PE) reads from each library.
	* **(5.a.2)** Mapping of the libraries to the indexed prokaryotic bins.
	* **(5.a.3)** Counting of the number of hits. 
	* **(5.a.4)** Generation of a crosstable with libraries as columns and prokaryotic bin representatives as rows and calculation of the unmapped reads percentage.
	* **(5.a.5)** Calculation of genome coverage using the BRAT results, the average number of base pairs per library and the number of base pairs and the total contigs per bin (optional).
* **5.b** Calculation of relative abundance and genome coverage of Uncultivated Viral Genomes and construction of relative abundance and genome coverage tables
	* **(5.b.1)** Merging of Paired-End (PE) reads from each library.
	* **(5.b.2)** Mapping of the libraries to the indexed viral contigs.
	* **(5.b.3)** Counting of the number of hits. 
	* **(5.b.4)** Generation of a crosstable with libraries as columns and viral contig representatives as rows and calculation of the unmapped reads percentage.
	* **(5.b.5)** Calculation of genome coverage using the BRAT results, the average number of base pairs per library and the number of base pairs and the total contigs per bin (optional).
* **5.c** Calculation of relative abundance  and genome coverage of Eukaryotic Metagenome-Assembled Genomes and construction of relative abundance and genome coverage tables
	* **(5.c.1)** Merging of Paired-End (PE) reads from each library.
	* **(5.c.2)** Mapping of the libraries to the indexed eukaryotic bins.
	* **(5.c.3)** Counting of the number of hits. 
	* **(5.c.4)** Generation of a crosstable with libraries as columns and eukaryotic bin representatives as rows and calculation of the unmapped reads percentage.
	* **(5.c.5)** Calculation of genome coverage using the BRAT results, the average number of base pairs per library and the number of base pairs and the total contigs per bin(optional).
* **5.d** Construction of combined abundance tables
	* **(5.d.1)** Construction of combined Relative Abundance Table.
	* **(5.d.2)** Construction of combined Coverage Abundance Table (optional).



# System requirements

**Warning:** complex dependencies ahead!

Hopefully, you won't have to worry too much about it. We designed an installation script that should take care of every dependency for you.
The software used during the pipeline are detailed described here: [Dependencies description](https://github.com/JotaKas/MuDoGeR/blob/master/dependencies_description.md).

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






# Simplified usage of the MuDoGeR

A more detailed tutorial for the MuDoGeR can be found in [Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md) file. In this file, instructions and examples can be found.

After all dependencies are in position then the MuDoGeR is ready to run. The run of **MuDoGeR** is easy and there is a main script that can wrap around all the **MuDoGeR**'s individual module. The individual modules can be called independently.

```
mudoger -h
	Usage: mudoger [module] --help
	Options:
	
	read_qc		Raw read QC 
	resources	Calculation of memory required by the MetaSpades, for assembling the good quality reads 
	assembly	Assembly 
	prokaryotic	Recovery of Prokaryotic Metagenome-Assembled Genomes 
	viral		Recovery of Uncultivated Viral Genomes
	eukaryotic 	Recovery of Eukaryotic Metagenome-Assembled Genomes
	relative_abund	Relative abundance calculation
		
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
	--g STR		host folder directory
	-c INT		minimum coverage (default=70)
	-i INT		minimum identity (default=95)
	
	-virfinder	Recovery of viral data with VirFinder 
	-virsorter	Recovery of viral data with VirSorter 	
	-vibrant	Recovery of viral data with VIBRANT
	-dereplication	Removal of replicate sequences
	-checkv		Quality control of dereplicated contigs with CheckV		
	-vcontact2	Taxonomy of dereplicated contigs with vContact2
	--wish		Host identification of dereplicated contigs with WIsH
	-mags_select	Selection of viral MAGs
```

# Citing

# Acknowledgements

	



 


	

