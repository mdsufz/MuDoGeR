 # Multi-Domain Genome Recovery (MuDoGeR)
 
![](https://github.com/mdsufz/MuDoGeR/blob/master/MUDOGER_FRAMEWORK.JPG)


The **Multi-Domain Genome Recovery (MuDoGeR)** framework (**Figure 1**) is a tool developed to help users to recover Metagenome-Assembled Genomes (MAGs as defined by Parks et al. (2017)) and Uncultivated Viral Genomes (UViGs as defined by  Roux (2019)) from dozens to hundreds of samples simultaneously. The **MuDoGeR** framework combines pipelines from different tools. The **MuDoGeR** starts with Module 1: **Pre-Processing**, which is divided in 3 steps: **1.a** **Raw Read Quality Control** and **1.b** **Resources calculation**, which feed the **1.c** **Assembly**. These assemblies will be used in the other modules. After pre-processing of the data, **MuDoGeR** is divided in 3 different branches: Module 2: **Recovery of Prokaryotic Metagenome-Assembled Genomes** (using **metaWrap**), Module 3: **Recovery of Uncultivated Viral Genomes** (using **VirSorter**, **VirFinder** and **VIBRANT**) and Module 4: **Recovery of Eukaryotic Metagenome-Assembled Genomes** (using **EukRep**). Furthermore, in **Module 5**: **Relative Abundance**, a strategy was developed for mapping the relative abundance of the selected MAGs in each library. Also, a step was added for bin_refinement of the selected Metagenome-Assembled Genomes from Prokaryotes, using **U-bin**. References of the used tools can be found in the end of the page.
 

## MuDoGeR Overview

* Instructions for using the **MuDoGeR** can be found in the following hyperlink: ![Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md).
* Detailed description of the **MuDoGeR** steps can be found in the following hyperlink: ![Module description](https://github.com/mdsufz/MuDoGeR/blob/master/module_description.md).
* Information about the system requirements of the **MuDoGeR** can be found in the following hyperlink: ![System requirements](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#system-requirements).
* Detailed instructions for the installation of the **MuDoGeR** tools can be found in the following hyperlink: ![Installation](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#installation).
* The simplified usage of the **MuDoGeR** can be found in the following hyperlink: ![Simplified usage of the MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#simplified-usage-of-the-mudoger).

### Module 1: Pre-Processing 

![](https://github.com/mdsufz/MuDoGeR/blob/master/Module1.JPG)

 The steps of Module 1 are shown in Figure 2 and they are excecuted with the scripts find in the following hyperlink: ![Pre-Processing](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-1-pre-processing).
 
 The steps of Module 1 can be summarized as following:

* **1.a**: Raw Read Quality Control.
* **1.b**: Calculation of resources.
    * **(1.b.1)** The k-mer of the good quality reads produced in **1.a** is calculated. The k-mer sizes that will be investigated are 33 and 55. The results of both k-mer 33 and 55 are combined in a single file.  
    * **(1.b.2)** The calculated k-mer is added to an equation that is used to estimate the amount of memory that **metaSPades** utilizes to assemble the good-quality reads.
* **1.c**: Assembly.

### Module 2: Recovery of Prokaryotic Metagenome-Assembled Genomes

![](https://github.com/mdsufz/MuDoGeR/blob/master/Module2.JPG)

The different steps of the Module 2 are shown in Figure 3 and they are excecuted with the scripts find in the following hyperlink: ![Pipeline for recovery of Prokaryotic Metagenome-Assembled Genomes](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-2-recovery-of-prokaryotic-metagenome-assembled-genomes).

 The steps of Module 2 can be summarized as following:
 
* **2.a**: Binning of Prokaryotic Metagenome-Assembled Genomes, bin_refinement, taxonomic classification,quality estimation and annotation of Prokaryotic bins
     * **(2.a.1)** Binning with **MaxBin2**, **metaBAT2**, and **CONCOCT**. 
     * **(2.a.2)** Dereplication of bins for prior bin_refinement. The completeness/contamination parameters have been set to 50%/10% for Bacteria and 40%/30% for Archaea. 
     * **(2.a.3)** Taxonomic classification of the genomic bins produced in **(2.a.2)** using **GTDB-Tk**.
     * **(2.a.4)** Generation of quality matrix of genomic bins produced in **(2.a.2)** using **CheckM**. 
     * **(2.a.5)** Filtering of genomic bins produced in **(2.a.4)**, by bin quality. The minimum quality for the filtering is set by default at 50. In this step, the user can also change the required quality (optional step).      
     * **(2.a.6)** Annotation of genomic bins produced in **(2.a.2)** with **PROKKA**.
* **2.b**: Selection of Prokaryotic Metagenome-Assembled Genomes Representatives.
	* **(2.b.1)** Grouping of Prokaryotic Metagenome-Assembled Genomes by taxonomy, using the .tsv file of **GTDB-Tk** results.
	* **(2.b.2)** Separation of the cluster produced in **(2.b.1)** by species, using ANI (Average Nucleotide Identity) splitter, with default ANI_distance 0.95.
	* **(2.b.3)** Selection of representative Metagenome-Assembled Genomes, from the clusters produced in **(2.b.2)**.
* **2.c**: Refinement of the bins produced in binning or/and the bins of the selected Representative Metagenome-Assembled Genomes using **U-bin** (optional step).
    
    
### Module 3: Recovery of Uncultivated Viral Genomes 

![](https://github.com/mdsufz/MuDoGeR/blob/master/Module3.JPG)

The steps of the Module 3  are shown in Figure 4 and they are excecuted in the scripts find in the following hyperlink: ![Pipelines for viral genomes recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-3-recovery-of-uncultivated-viral-genomes).

 The steps of Module 3 can be summarized as following: 
 
* **3.a**: Recovery, quality estimation, taxonomic classification and host identification of Uncultivated Viral Genomes 
    * **(3.a.1)** Recovery of Uncultivated Viral Genomes using **VirFinder**, **VirSorter** and **VIBRANT**.
    * **(3.a.2)** Filtering of the recovered genomes. From the results of **VirFinder**, the sequences with p-value > 0.01 and length < 1000 bp are removed. From the results of **VirSorter**, which recovers prophages, only the sequences of categories 1 and 2 are kept. From the results of **VIBRANT**, the combined assemblies of the phages are kept.
    * **(3.a.3)** Combination of the headers of all the Uncultivated Viral Genomes in a single file. Then, removal of the headers of the repeated sequences and sorting the remaining headers by sequence length.
    * **(3.a.4)** Assembly of the headers in the file produced in **1.c** with the respective sequences from the assembly file produced by **metaSPades** in **1.c**, for the generation of a fasta file.
    * **(3.a.5)** Removal of replicates from the assemblies in the extracted file using **Stampede-clustergenomes** with minimum coverage of 70% and minimum identity of 95%.
    * **(3.a.6)** Checking the quality of the dereplicated contigs with **CheckV**.
    * **(3.a.7)** Taxonomic classification of the clean contigs produced by **CheckV**,  with **vContact2**. In small datasets, we can use the dereplicated contigs from step **(3.a.5)**.
    * **(3.a.8)** Host identification of the clean contigs produced by **CheckV**, using **WIsH**. In small datasets, we can use the dereplicated contigs from **(3.a.5)** (optional).    
* **3.b**: Selection of Uncultivated Viral Genomes 
	* **(3.b.1)** Selection of all viruses that yielded taxonomy when using vContact2 plus those larger than 15 Kb.
   

### Module 4: Recovery of Eukaryotic Metagenome-Assembled Genomes

![](https://github.com/mdsufz/MuDoGeR/blob/master/Module4.JPG)

The steps of the Module 4  are shown in Figure 5 and they are excecuted with the scripts find in the following hyperlink: ![Pipelines for eukaryotic genomes recovery](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-4-recovery-of-eukaryotic-metagenome-assembled-genomes).

 The steps of Module 4 can be summarized as following:
* **4.a**: Recovery of Eukaryotic assemblies and production of Eukaryotic bins.
    * **(4.a.1)** Classification of Eukaryotic assemblies and removal of prokaryotic assemblies with **EukRep**.
    * **(4.a.2)** Use of **CONCOCT** for binning the Eukaryotic assemblies.  
    * **(4.a.3)** Filtering the Eukaryotic bins, produced from **CONCOCT**, by size. Bins with size < 2.5 Mb are removed.
* **4.b**: Completeness/contamination estimation and annotation of Eukaryotic bins
    * **(4.b.1)** In the filtered bins produced in **4.a**, genes are predicted by the trained **GeneMark-ES** model.
    * **(4.b.2)** **EukCC** utilization for estimating the contamination of Eukaryotic filtered bins produced in **4.a**.
    * **(4.b.3)** **MAKER2** annotates the predicted genes produced by **GeneMark-ES**. 
    * **(4.b.4)** **BUSCO** is applied to the annotated genes from **MAKER2**, for detection of single copy orthologous genes (SCGs) and estimation of the completeness of Eukaryotic contigs.
* **4.c**: Selection of Eukaryotic Metagenome-Assembled Genomes Representatives.
	* **(4.c.1)** Grouping of Eukaryotic Metagenome-Assembled Genomes by taxonomy, using **BUSCO** and **EukCC** results.
	* **(4.c.2)** Separation of the cluster produced in **(4.c.1)** by species, using ANI (Average Nucleotide Identity) splitter, with deafault ANI_distance 0.95.
	* **(4.c.3)** Selection of representative Metagenome-Assembled Genomes, from the clusters produced in **(4.c.2)**.

### Module 5 Relative abundace
![](https://github.com/mdsufz/MuDoGeR/blob/master/module5.JPG)

The steps of the Module 5 are shown in Figure 6 and they are excecuted with the scripts find in the following hyperlink: ![Pipelines for Relative abundace ](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#module-5-relative-abundance).
* **5.a** Calculation of relative abundance and genome coverage of Prokaryotic Metagenome-Assembled Genomes and construction of relative abundance and genome coverage tables
	* **(5.a.1)** Merging of Paired-End (PE) reads from each library.
	* **(5.a.2)** Mapping of the libraries to the indexed prokaryotic bins.
	* **(5.a.3)** Counting of the number of hits. 
	* **(5.a.4)** Generation of a crosstable with the libraries in the columns and the prokaryotic bins in the rows and calculation of the unmapped reads percentage.
	* **(5.a.5)** Calculation of genome coverage using the BRAT results, the average number of base pairs per library and the number of base pairs and the total contigs per bin.
* **5.b** Calculation of relative abundance and genome coverage of Uncultivated Viral Genomes and construction of relative abundance and genome coverage tables
	* **(5.b.1)** Merging of Paired-End (PE) reads from each library.
	* **(5.b.2)** Mapping of the libraries to the indexed viral contigs.
	* **(5.b.3)** Counting of the number of hits. 
	* **(5.b.4)** Generation of a crosstable with the libraries in the columns and the viral contigs in the rows and calculation of the unmapped reads percentage.
	* **(5.b.5)** Calculation of genome coverage using the BRAT results, the average number of base pairs per library and the number of base pairs and the total contigs per bin.
* **5.c** Calculation of relative abundance  and genome coverage of Eukaryotic Metagenome-Assembled Genomes and construction of relative abundance and genome coverage tables
	* **(5.c.1)** Merging of Paired-End (PE) reads from each library.
	* **(5.c.2)** Mapping of the libraries to the indexed eukaryotic bins.
	* **(5.c.3)** Counting of the number of hits. 
	* **(5.c.4)** Generation of a crosstable with the libraries in the columns and the eukaryotic bins in the rows and calculation of the unmapped reads percentage.
	* **(5.c.5)** Calculation of genome coverage using the BRAT results, the average number of base pairs per library and the number of base pairs and the total contigs per bin.
* **5.d** Construction of combined abundance tables
	* **(5.d.1)** Construction of combined Relative Abundance Table.
	* **(5.d.2)** Construction of coverage Abundance Table (optional).


	

# System requirements
MuDoGeR is designed to support only Linux x64 systems. As for the resource requirements,  MuDoGeR framework uses many software that require large amount of memory (e.g **GDTB-Tk**, **MetaWRAP** ). Also, the resource requirements varies depending on the amount of the used data. Based on that, we recommend the user to use at least 180 GB RAM. The user also should take into account the amount of memory required for the assembly. As mentioned above, in the **1.c** step, the user can calculate the amount of memory for **metaSPades** assembly. The user should be aware that samples with higher diversity require higher amount of memory. In the following session, the dependencies for the tools used in the MuDoGeR framework are listed. 

## Dependencies
### Metawrap

* **Metawrap** uses python version of 2.7.14.

| Dependencies |  
|---|
| Database Checkm_DB |   
| Database KRAKEN standard database |   
| Database KRAKEN2 standard database |  
| Database NCBI_nt|   
| Database NCBI_tax | 
| Database Indexed hg38| 
| GCC/6.4.0-2.28 | 
| OpenMPI/2.1.2 | 

The Databases can be installed by following the instructons of the following hyperlink: ![Databases_metaWrap](https://github.com/bxlab/metaWRAP/blob/master/installation/database_installation.md).

### GTDB-Tk

* **GTDB-Tk** uses Python version >=3.6 (3.6.4).


| Dependencies | Version |
|---|---|
| Library OpenMPI | 2.1.2 |
| Package GCC | 6.4.0-2.28 |
| Lirbrary DendroPy | >= 4.1.0 |
| Lirbrary NumPy | >= 1.9.0 |
| Lirbrary tqdm | >= 4.31.0  |
| Software Prodigal | >= 2.6.2 |
| Software HMMER | >= 3.1b2 |
| Software pplacer | >= 1.1 |
| Software FastANI | >= 1.32 |
| Software FastTree | >= 2.1.9 |
| Software Mash |>= 2.2 |

* The dependencies of **GTDB-Tk** can be found in the following webpage: https://ecogenomics.github.io/GTDBTk/installing/index.html.

### CheckM

* **CheckM** uses Python version 3.7.4.

| Dependencies | Version |
|---|---|
| Package GCCcore | 8.3.0 |
| Software HMMER | >=3.1b1 |
| Software Prodigal | 2.60 or >=2.6.1 |
| Software pplacer | >=1.1 |
| Lirbrary NumPy | >= 1.9.0 |
| Lirbrary tqdm | >= 4.31.0  |
| Software Prodigal | >= 2.6.2 |

* As pplacer requires also dependencies, these dependencies, together with instructions of how to install them, can be found in the following link: http://matsen.github.io/pplacer/compiling.html.


### Ubin (optional)
* The pre-requests for **uBin** and how to install them, can be found in the follwoing hyperlink: ![Creating input files for uBin](https://github.com/ProbstLab/uBin-helperscripts#creating-input-files-for-ubin).

### DAS_Tool (optional)
* The dependencies of **DAS_Tool** can be found in the following hyperlink: ![Dependencies_DAS_TOOL](https://github.com/cmks/DAS_Tool#dependencies).

### VirSorter 

| Dependencies | 
|---|
| Software HMMER of version 3.1b2 | 
| Software MCL |  
| Software Metagene Annotator |
| Software MUSCLE | 
| Software BLAST+ | 
| Software DIAMOND |
| Package gcc/4/8.1-3 | 

* The software used in the **VirSorter** can be found also in the following hyperlink: ![Software_VirSorter](https://github.com/simroux/VirSorter#dependencies).

### CheckV

* For **CheckV** to work the installation of Anaconda3/5.3.0 with python version of 3.6 is required.

| Dependencies | Version  |  
|---|---|
| Software DIAMOND | 2.0.4 | 
| Software HMMER | 3.3 | 
| Software Prodigal | 2.6.3 |  

* More information about **CheckV** dependencies can be found by following the webpage: https://bitbucket.org/berkeleylab/checkv/src/master/.

### vContact2

* **vContact2** uses Anaconda3/5.3.0 and python of version >=3.7 (3.7.4 ).

| Dependencies | Version | 
|---|---|
| Package GCCcore | 8.3.0 |
| Package Java/11.0.2 |
| biopython | >=1.73 |  
| hdf5 |>=1.10.4 |  
| networkx | >=2.2 | 
| numpy | >=1.15.4 | 
| pandas | >=0.25.0,<=0.25.3 | 
| psutil | >=5.5.0 | 
| pyparsing  | >=2.4.6 |
| pytables |  >=3.4.0 | 
| scikit-learn | >=0.20.2 | 
| scipy | >=1.2.0 |

*  More information about **vContact2** dependencies can be found by following the webpage: https://bioconda.github.io/recipes/vcontact2/README.html.

### VIBRANT

* **VIBRANT** uses python of version >= 3.5 (3.6.6).

| Dependencies | 
|---|
| KEGG | 
| Pfam (v32) | 
| VOG (release 94) |  
| Package GCC/7.3.0-2.30 |  
| Package OpenMPI/3.1.1 |  
| Prodigal version 2.6.3 |    
| HMMER3 |  
| gzip | 
| tar | 
| wget |
| BioPython | 
| Pandas |  
| Matplotlib   |
| Seaborn |   
| Numpy (version >= 1.17.0) |  
| Scikit-learn (version == 0.21.3)|
| Pickle |

* More information about **VIBRANT** dependencies can be found in the following hyperlink: ![VIBRANT_Dependencies](https://github.com/AnantharamanLab/VIBRANT#requirements-).

### Stampede-clustergenomes
* **Stampede-clustergenomes** uses python of version 2.7.14.

| Packages | Version |    
|---|---| 
| MUMmer |  3.0 |
| GCCcore | 6.4.0 |  
| libtool | 2.4.6 |  
| Perl | 5.26.0 |  
| OpenMPI |2.1.2 | 
| GLib | 2.54.3 | 

* Instructions for installation of Perl/5.26.0 and MUMmer/3.0 can be found in the following webpage: https://bitbucket.org/MAVERICLab/stampede-clustergenomes/src/master/.

### VirFinder

* **VirFinder** uses R/3.4.4. Also the user need to have OpenMPI/2.1.2 library installed.    

| Dependencies | 
|---|
| R package glmnet |  
| R package Rcpp |  
| R package qvalue |   
| Package GCC/6.4.0-2.28 | 
| Package OpenMPI/2.1.2 | 

### WIsH

* **WIsH** requires OpenMP library installed.

### EukRep

* **EukRep** uses python of version 3.6.4.
 
| Dependencies | 
|---|
| GCC/6.4.0-2.28 |

### MAKER2

| Dependencies |  
|---|
| RepeatModeler(1.0.4) |
| RepeatMasker (4.0.5) |
| GCC/6.4.0-2.28 |
| RepBase (version used was 20150807) | 
| MAKER/2.31.10-foss-2019b-1 |  
| Augustus version 3.3 |  
| BUSCO version 3 |  
| SNAP |
| BEDtools version 2.24.0 |

* As SNAP requires some more dependencies, they can be found in the following hyperlink: ![Compiling_snap](https://github.com/KorfLab/SNAP#compiling).
* As RepeatModeler(1.0.4)  requires some more dependencies,they can be found in the following hyperlink: ![RepeatModeler_dependencies](https://github.com/Dfam-consortium/TETools).
* As RepeatMasker (4.0.5) requires some more dependencies, they can be found in the following link: http://www.repeatmasker.org/RepeatMasker/.
* As Augustus requires some more dependencies, they can be found in the following hyperlink: ![Dependencies_Augustus](https://github.com/Gaius-Augustus/Augustus#install-dependencies).
* As BEDtools requires some more dependencies, they can be found in the following link: https://bedtools.readthedocs.io/en/latest/content/installation.html.
* About Busco, the user can check the following hyperlink:  ![BUSCO_Manual](https://github.com/WenchaoLin/BUSCO-Mod).
* More information about **MAKER2** dependencies can be found in the following hyperlink: ![Maker2_Software prerequisites](https://github.com/wuying1984/MAKER2_PM_genome_annotation#software-prerequisites).

### BUSCO
* **BUSCO** uses python of version 3.7.4.

| Dependencies |  
|---|
| Package GCC/6.4.0-2.28 | 
| Package OpenMPI/2.1.2 |  
| foss/2019b |  
| NCBI BLAST+ |
| HMMER 3.1b2 |
| Augustus 3.0.x |

* As Augustus requires some more dependencies, they can be found in the following hyperlink: ![Dependencies_Augustus](https://github.com/Gaius-Augustus/Augustus#install-dependencies).
* More information about **BUSCO** dependencies can be found in the following hyperlink:![BUSCO_Dependencies](https://github.com/WenchaoLin/BUSCO-Mod#setup).

### EukCC/GeneMark-ES
* The **EukCC** uses python of version 3.7.4. Also the **EukCC** needs also installation of the **GeneMark-ES** first which has more dependencies. For the installation of **GeneMark-ES** dependencies the user can run the following command which will install all of them:

```
apt install -y cpanminus make gcc  dialog

cpanm inc::Module::Install::DSL Hash::Merge MCE::Mutex FindBin Test::Pod Logger::Simple  Parallel::ForkManager.pm YAML
```

| Dependencies for **EukCC/GeneMark-ES** |  
|---|
| Package GCCcore/8.3.0 | 
| Package foss/2019b  |  
| Perl/5.30.0 |
| Anaconda2/5.3.0 |
| Package foss-2019b |

* More information about **EukCC/GeneMark-ES** dependencies can be found by following the webpage: https://eukcc.readthedocs.io/en/latest/install.html#install-via-conda.

### Ani_Splitter

| Dependencies |  
|---|
| Package GCCcore/8.3.0 | 
| Package OpenMPI/3.1.4 |  
| R |
| cluster | 
| dendextend |  
| dplyr |
| fpc | 
| ggplot2 |  
| reshape2 |
| rlang | 
| tidyselect |  
| vctrs |



### BRAT 

| Dependencies | 
|---|
| Python 3.6.6 |
| OpenMPI 3.1.1 |
| GCC 7.3.0-2.30 |
| datamash-1.3 |
| Bowtie2 2.3.5.1 |
| SAM tools 1.9 |

* More information about Bowtie 2 can be found by following the webpage: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml.
* More information about Samtools can be found by following the webpage: http://www.htslib.org/doc/samtools.html.
* More information about can be found by following the webpage: http://brat.nlplab.org/installation.html.

### For resource calculation

| Dependencies |  
|---|
| Package GCC/8.3.0 |
| Package OpenMPI/3.1.4 |
| R/3.6.2-2 |






# Installation

## Pre-Processing Module
### TrimGalore
The instructions for  TrimGalore installation can be found at ![TrimGalore_Manual](https://github.com/FelixKrueger/TrimGalore).
### Bmtagger
The instructions for the bmtagger installation can be found by following the webpage: https://bioconda.github.io/recipes/bmtagger/README.html.
### MegaHit
The instructions for the MegaHit installation can be found at ![MegaHit_Manual](https://github.com/voutcn/megahit).
### MetaSPades
The instructions for the metaSPAdes installation can be found at ![metaSPAdes_Manual](https://github.com/ablab/spades).
### FastQC
Instructions for FastQC installation can be found at ![FastQC_Manual](https://github.com/pnnl/fqc).

## Prokaryotic module
### MetaWRAP
The instructions for the metaWrap installation can be found at ![metaWRAP_Manual](https://github.com/bxlab/metaWRAP).
### MaxBin2
The instructions for the MaxBin2 installation can be found at ![MaxBin2_Manual](https://github.com/assemblerflow/flowcraft/blob/master/docs/user/components/maxbin2.rst).
### MetaBAT2
The instructions for the  metaBAT2 installation can be found by following the webpage: https://bioconda.github.io/recipes/metabat2/README.html.
### CONCOCT 
The instructions for the CONCOCT installation can be found at ![CONCOCT_Manual](https://github.com/BinPro/CONCOCT).
### CheckM
The instructions for the CheckM  installation can be found at ![CheckM_Manual](https://github.com/Ecogenomics/CheckM).
### uBin
The instructions for the uBin installation can be found at ![uBin_Manual](https://github.com/ProbstLab/uBin).
### GTDB-Tk
The instructions for the GTDB-Tk installation can be found at ![GTDB-Tk_Manual](https://github.com/ecogenomics/gtdbtk).
### DAS_Tool
The instructions for the DAS_Tool installation can be found at ![DAS_Tool_Manual](https://github.com/cmks/DAS_Tool).

## Viral module
### VirFinder
The instructions for VirFinder installation can be found at ![VirFinder_Manual](https://github.com/jessieren/VirFinder).
### VirSorter
The instructions for VirSorter installation can be found at ![VirSorter_Manual](https://github.com/simroux/VirSorter). 
### VIBRANT 
The instructions for VIBRANT installation can be found at ![VIBRANT_Manual](https://github.com/AnantharamanLab/VIBRANT).
### Stampede-clustergenomes
The instructions for Stampede-clustergenomes usage cane be found by following the webpage: https://bitbucket.org/MAVERICLab/stampede-clustergenomes/src/master.
### CheckV
The instructions for CheckV installation can be found by following the webpage: https://bitbucket.org/berkeleylab/checkv/src/master/.
### vContact2
The instructions for vContact2 installation can be found by following the webpage: https://bioconda.github.io/recipes/vcontact2/README.html.
### WIsH
The instructions for WIsH installation ![WIsH_Manual](https://github.com/soedinglab/wish).

## Eukaryotic module
### EukRep
The instructions for EukRep installation can be found at ![EukRep_Manual](https://github.com/patrickwest/EukRep_Pipeline ). 
### EukCC/GeneMark
The instructions for EukCC/GeneMark installation can be found at ![EukCC/GeneMark_Manual](https://github.com/Finn-Lab/EukCC).
### BUSCO
The instructions for BUSCO installation can be found at ![BUSCO_Manual](https://github.com/WenchaoLin/BUSCO-Mod).
### MAKER2
The instructions for MAKER2 installation can be found at ![MAKER2_Manual](https://github.com/wuying1984/MAKER2_PM_genome_annotation).

## Relative abundance
### BRAT
The instructions for BRAT installation can be found at  ![BRAT_Manual](https://github.com/nlplab/brat).

# Simplified usage of the MuDoGeR

A more detailed tutorial for the MuDoGeR can be found in ![Manual_MuDoGeR](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md) file. In this file, instructions and examples can be found.

After all dependencies are in position then the MuDoGeR is ready to run. The run of **MuDoGeR** is easy and there is a main script that can wrap around all the **MuDoGeR**'s individual module. The individual modules can be called independently.

```
mudoger -h
	Usage: mudoger [module] --help
	Options:
	
	read_qc		Raw read QC 
	resources	Calculation of memory required by the MetaSpades, for assembling the good quality reads 
	assembly	Assembly 
	prokaryotic	Recovery of Prokaryotic etagenome-Assembled Genomes 
	viral		Recovery of Viral Metagenome-Assembled Genomes
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

	



 


	

