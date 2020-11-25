# DETAILED DESCRIPTION OF MuDoGeR

## Module 1: Pre-Processing 
Module 1 is devided into 3 steps:

### 1.a: Raw Read Quality Control
The **MuDoGeR** Raw Read Quality Control is based on **metaWrap** Read Quality Control and is responsible to prepare the raw reads for assembly and alignment. The trimming of the raw reads is based on adapted content and **PHRED** scored with the default setting of Trim-galore. As a result, all the non-high quality sequences are removed. The high quality reads that remained are aligned to the genome of the host with bmtagger, a process that leads to the removal of host contamination and of read pairs with only one aligned read from the metagenomic data. In the final step, **FASTQC** is used to assess read quality improvement. Full detailed description of the Read_qc can be found in the following hyperlink: ![metaWrap Read_qc](https://github.com/bxlab/metaWRAP/blob/master/Module_descriptions.md#read_qc).

### 1.b: Calculation of resources 
Before the user goes to the assembly module, it is possible to estimate the necessary resources for the assembly of the library produced in quality control module. This can be achieved by estimating the complexity of the forward or reveresed pure reads by k-mer calculation, as they are almost identical. The k-mer is usually estimated in size 33 or 55. Both values should be estimated and  k-mer with sizes of 33 and 55 are combined to a single output file **(1.b.1)**. Using this file the necessary amount of memory that **metaSPades** needs to utilze for assembly the good qulity reads is predicted **(1.b.2)**.

### 1.c: Assembly 
The **MuDoGeR** Assembly is based on **metaWrap** Assembly and is utillized for assembling a set of metagenomic reads. Two possible readers can be chosen for the assembly: **MegaHit** (default) and **metaSPAdes**. **MegaHit** reader requires less memory and is also faster, compared to **metaSPAdes**. On the other hand, **metaSPAdes** reader has been proved more reliable as it produces assemblies of higher quality. Because of that, the usage of **metaSPAdes** reader is preferable except of the cases of large data-set. Fully detailed description of the metaWrap Assembly can be found in the following hyperlink: ![metaWrap Assembly](https://github.com/bxlab/metaWRAP/blob/master/Module_descriptions.md#assembly).  

## Module 2 Recovery of Prokaryotic Metagenome-Assembled Genomes
Module 2 is separated in tasks:

### 2.a: Binning of Prokaryotic Metagenome-Assembled Genomes, bin_refinement, quality estimation, taxonomical classification and annotation of Prokaryotic bins

This step is a combination of tasks included in the **metaWrap** tool. The module starts with the binning of the assembly data-sets with the use of three wrapped tools: **MaxBin2**, **metaBAT2** and **CONCOCT** **(2.a.1)**. This process results to the production of bin files from each method. Next, the bin_refinement of the bins produced from each tool follows **(2.a.2)**. In the bin_refinement, the minimum completeness and contamination for bacterial bins are 50% and 10% respectively, while the minimum completeness and contamination for the archael bins are 40% and 30% respectively. For the taxonomic classification of the bin sets from Bacteria and Archaea produced in the bin refinement module, **GTDB-Tk** tool is used **(2.a.3)**. Furthermore, the quality matrix of the bins produced from bin_refinement, using **CheckM** **(2.a.4)**. After this step, the user can choose to filter the **CheckM** bins by quality **(2.a.5)**. Bin quality is defined as completeness – 5×contamination (Parks, 2018). By default, the quality value is 50, but the user can choose a different value. Finally, the prediction of functional annotation of the bins produced in the bin_refinement is achieved by **PROKKA**, which utilizes a number of softawre tools **(2.a.6)**. Full detailed description of **(2.a)** can be found in the following hyperlink: ![Detailed descriptions of each module](https://github.com/bxlab/metaWRAP/blob/master/Module_descriptions.md). 

### 2.b: Selection of Prokaryotic Metagenome-Assembled Genomes, uBin-refinement and relative abundance (not done, so maybe more will written)

In this step, a selection of Prokaryotic Metagenome-Assembled Genomes takes place **(2.b.1)**. In case the user chooses to, **DAS Tool** and **uBin** refining tool are used for the manual curation of genomes from the selected metagenomic bins **(2.b.1)**. Finally the relative abundance table is contructed **(2.b.3)**. 

## Module 3: Recovery of Viral Metagenome-Assembled Genomes
Module 3 is devided into 2 steps:

### **3.a**: Recovery, quality estimation, taxonomical classification and host identification of Viral Metagenome-Assembled Genomes

In the beginning of this step, three viral recovery tools (**VirFinder**, **VirFinder**, **VIBRANT**) are utilized for the identification and recovery of viral genomes from an assemble metagenomic data-set. In the beginning, each of the tools is used to recover independently the viral genomes from a given assembled data-set and saves them in seperate folders **(3.a.1)**.  Following that, the recovered sequences are filtered **(3.a.2)**. The selection of proper **VirFinder** sequences is based on low q-value (q-value=< 0.01) and high length (length>= 1000 bp), the **VirSorter** chosen sequences are those classified into categories 1 and 2 and from **VIBRANT** the selected sequences are those of the combined assemblies from phages. The headers of **VirSorter** and **VIBRANT** contigs are modified so they can match with those produced from the **VirFinder** results. It is important to note that in contrast with the other two tools, the **VirFinder** output recovery file contains only the headers of the assemblies. Because of that, the headers of the unique filtered sequences from each tool are extracted to a common fasta file and sorted by length **(3.a.3)**. Then, using the headers including in the common fasta file, the actual sequences from the assembly data-set are extracted and transfered to a new fasta file **(3.a.4)**. Next,the duplicated contigs are removed by dereplication using **Stampede-clustergenomes** tool **(3.a.5)**. The dependencies of the dereplication step are the maximum coverage (-c) 70% and the maximum identity (-i) 95%. The dereplication is followed by quality estimation (completeness – 5×contamination (Parks, 2018)) **(3.a.6)** and taxonomy **(3.a.7)** of the the dereplicated contigs by **CheckV** and **vContact2** respectively. Finally, if the user chooses to, the hosts of the dereplicated contigs are identified by using **WIsH tool** **(3.a.8)**.

### **3.b**: Selection of Viral Metagenome-Assembled Genomes, uBin-refinement and relative abundance (not done, so maybe more will written)

In this step a selection of Prokaryotic Metagenome-Assembled Genomes takes place **(3.b.1)** and the relative abundance table is constructed **(3.b.2)**.

## Module 4: Recovery of Eukaryotic Metagenome-Assembled Genomes
Module 4 is devided into 3 steps:

### 4.a: Recovery of Eukaryotic assemblies and production of Eukaryotic bins
In the beginning this step, the assembled data-set is separated to Prokaryotic and Eukaryotic assemblies using **EukRep** tool. In case the size of the Eukaryotic assembly file is >= 2.5 Mb, the Eukaryotic recovery modules can continue to the automated binnning with **CONCOCT** tool. The eukaryotic bins  produced by **CONCOCT**. Finally. the bins are filtered by size and those of size < 2.5 Mb are removed. 

### 4.b: Completeness/contamination estimation and annotation of Eukaryotic bins

In this step, a chain of processes is followed for one of the bins produced in **4.b**: In the beginning, **GeneMark-EV** is applied for gene prediction **(4.b.1)**. Also, the contamination of the bins which were kept after the filtering is estimated by using **EukCC** tool **(4.b.2)**. Then the predicted genes from **GeneMark-EV** are annotated with **Maker2** **(4.b.3)**. Finally the completeness of the annotated proteins is measured using **BUSCO** **(4.b.4)**. 

Finally, a selection of eukaryotic MAGs takes place.

