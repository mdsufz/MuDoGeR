# DETAILED DESCRIPTION OF MuDoGeR

## Module 1: Pre-Processing 
Module 1 is divided into 3 steps:

### 1.a: Raw Read Quality Control
The **MuDoGeR** Raw Read Quality Control is based on **metaWrap** Read Quality Control (Uritskiy et al., 2018) and is responsible to prepare the raw reads for assembly and alignment. The trimming of the raw reads is based on adapted content and **PHRED** scored with the default setting of Trim-galore (Krueger, 2015). As a result, all the non-high-quality sequences are removed. The high-quality reads that remained are aligned to the genome of a potential host with bmtagger (Agarwala and Morgulis, 2010), a process that leads to the removal of host contamination and of read pairs with only one aligned read from the metagenomic data. In the final step, **FASTQC** (Brown et al., 2017) is used to assess read quality improvement.

### 1.b: Calculation of resources
In case the user is interested in a fine quality assembly, then **metaSPades** (Nurk et al., 2017) should be used. Due to the high memory requirements of **metaSPades**, prior to the assembly, it is possible to estimate the necessary resources for the assembly of the library produced in **1.a**. **(1.b.1)** This can be achieved by estimating the complexity of the forward or reversed pure reads by k-mer calculation, as they are almost identical. The k-mer is usually estimated in size 33 or 55. Both k-mer values are estimated and k-mer with sizes of 33 and 55 are combined to a single output file. **(1.b.2)** Using this file, the necessary memory that **metaSPades** needs for the assembly of the good-quality reads is predicted.

### 1.c: Assembly
The **MuDoGeR** assembly is based on **metaWrap** assembly and is utillized for assembling a set of metagenomic reads. Two possible tools can be chosen for the assembly: **MegaHit** (default) and **metaSPAdes**. **MegaHit** (Li et al., 2016) requires less memory and is also faster, compared to **metaSPAdes** (if the user decides to utilize **metaSPAdes** we recommend the use of **1.b**). If the user is interested in a fine quality assembly, he will use **metaSPAdes** and not **MegaHit**. Due to the higher memory requirements of metaspades prior of the assembly it is possible to estimate the resources necessary for metaSPades usage. Because of that, the usage of **metaSPAdes** reader is preferable except of the cases of large data-set.

## Module 2: Recovery of Prokaryotic Metagenome-Assembled Genomes
Module 2 is divided in 3 steps:

### 2.a: Binning of Prokaryotic Metagenome-Assembled Genomes, bin_refinement, taxonomic classification, quality estimation and annotation of Prokaryotic bins
This step is a combination of tasks included in the **metaWrap** tool (Uritskiy et al., 2018). **(2.a.1)** This step starts with the binning of the assembly data-sets with the use of three wrapped tools: **MaxBin2** (Wu et al., 2016), **metaBAT2** (Kang et al., 2019) and **CONCOCT** (Alneberg et al., 2014). This process results in the production of bin files from each method. Next, **(2.a.2)** the bin_refinement of the bins produced from each tool follows. In the bin_refinement, the minimum completeness and contamination for bacterial bins are 50% and 10% respectively, while the minimum completeness and contamination for the archaeal bins are 40% and 30% respectively. **(2.a.3)** For the taxonomic classification of the bin sets from Bacteria and Archaea produced in the bin refinement module, **GTDB-Tk** tool (Chaumeil et al., 2020) is used. Furthermore, **(2.a.4)** the quality matrix of the bins produced from bin_refinement, using **CheckM** (Parks et al., 2015). After this step, **(2.a.5)** the user can choose to filter the **CheckM** bins by quality. Bin quality is defined as completeness – 5×contamination (Parks et al., 2018). By default, the quality value is 50, but the user can choose a different value. Finally, **(2.a.6)** the prediction of functional annotation of the bins produced in the bin_refinement is achieved by **PROKKA** (Seemann, 2014), which utilizes a number of software tools.

### 2.b: Selection of Prokaryotic Metagenome-Assembled Genomes Representatives
In this step, a selection of Prokaryotic Metagenome-Assembled Genomes according Nayfach et al. (2020) takes place. It uses average nucleotide identities (ANI) distances and the threshold is 0.95 (species-level OTUs). **(2.b.1)** First, the Prokaryotic Metagenome-Assembled Genomes are grouped by taxonomy. **(2.b.2)** Furthermore, the produced cluster are separated by species using ANI (Average Nucleotide Identity) splitter with default ANI_distance 0.95. Clusters with bootstrap > 75 are represent new species. **(2.b.3)** Finally, there is a selection of representative Metagenome-Assembled Genomes, from the clusters produced in **(2.b.2)**.

### 2.c: Refinement of the bins produced in binning or/and the bins of the selected Representative Metagenome-Assembled Genomes using U-bin
In case the user chooses to, **uBin** (Bornemann et al., 2020) refining tool can be used instead of **DAS Tool** (Sieber et al., 2018) in the bins produced in **(2.a.1)**. The users have also the option of using **uBin** to manually refined the bins produced in **2.b**. The hyperlink of the **uBin** tool can be found in the following hyperlink: ![ubin](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#ubin).

## Module 3: Recovery of Uncultivated Viral Genomes
Module 3 is divided into 2 steps: 

### **3.a**: Recovery, quality estimation, taxonomic classification and host identification of Uncultivated Viral Genomes
In the beginning of this step, three viral recovery tools (**VirFinder** (Ren et al., 2017), **VirSorter** which recovers prophages (Roux et al., 2015) and **VIBRANT** (Kieft et al., 2020)) are utilized for the identification and recovery of Uncultivated Viral genomes (Roux, 2019), from an assemble metagenomic data-set. **(3.a.1)** Each of the tools is used to recover independently the viral genomes from a given assembled data-set and saves them in separate folders. **(3.a.2)** Following that, the recovered sequences are filtered. The selection of proper **VirFinder** sequences is based on low q-value (q-value=< 0.01) and high length (length>= 1000 bp), the **VirSorter** chosen sequences are those classified into categories 1 and 2 and from **VIBRANT** the selected sequences are those of the combined assemblies from phages. **(3.a.3)** The headers of **VirSorter** and **VIBRANT** contigs are modified so they can match with those produced from the **VirFinder** results. It is important to note that in contrast with the other two tools, the **VirFinder** output recovery file contains only the headers of the assemblies. Because of that, the headers of the unique filtered sequences from each tool are extracted to a common fasta file and sorted by length. **(3.a.4)** Using the headers including in the common fasta file, the actual sequences from the assembly data-set are extracted and transferred to a new fasta file. **(3.a.5)** Next, the duplicated contigs are removed by dereplication using **Stampede-clustergenomes** tool (Bolduc and Roux, 2017). The dependencies of the dereplication step are the minimum coverage (-c) 70% and the minimum identity (-i) 95%. The dereplication is followed by **(3.a.6)**  quality estimation  of the dereplicated contigs by **CheckV** (Nayfach et al., 2020). **(3.a.7)** The taxonomic classification of the clean contigs produced by **CheckV** is achieved by **vContact2** (Jang et al., 2019). In case of small datasets, we can use the results of the dereplication process. If the user chooses to, **(3.a.8)** the hosts of the clean contigs, or in case of small datasets, the dereplicated contigs, are identified by using **WIsH tool** (Galiez et al., 2017). 

### **3.b**: Selection of Uncultivated Viral Genomes
In **3.b** the user can select viral representatives. **(3.b.1)** The viral representatives encompass all viruses that yielded taxonomy when using vContact2 plus those larger than 15 Kb.


## Module 4: Recovery of Eukaryotic Metagenome-Assembled Genomes
Module 4 is divided into 3 steps:

### 4.a: Recovery of Eukaryotic assemblies and production of Eukaryotic bins
In the beginning this step, **(4.a.1)** the assembled data-set is separated to Prokaryotic and Eukaryotic assemblies using **EukRep** tool (West et al., 2017). **(4.a.2)** In case the size of the Eukaryotic assembly file is >= 2.5 Mb, the Eukaryotic recovery modules can continue to the automated binning with **CONCOCT** tool. The eukaryotic bins  produced by **CONCOCT**. Finally, **(4.a.3)** the bins are filtered by size and those of size < 2.5 Mb are removed. 

### 4.b: Completeness/contamination estimation and annotation of Eukaryotic bins
In this step, a chain of processes is followed for one of the bins produced in **4.a**: **(4.b.1)** **GeneMark-ES** is applied for gene prediction. **(4.b.2)** Also, the contamination of the bins which were kept after the filtering is estimated by using **EukCC** tool (Saary et al., 2020). **(4.b.3)** The predicted genes from **GeneMark-ES** (Saary et al., 2020) are annotated with **MAKER2** (Holt et al., 2011). **(4.b.4)** The completeness of the annotated proteins is measured using **BUSCO** (Simão et al., 2015). 

### 4.c: Selection of Eukaryotic Metagenome-Assembled Genomes Representatives
Currently, a selection of Eukaryotic Metagenome-Assembled Genomes Representatives has not yet been benchmarked. Therefore, we are using standards benchmarkers with prokaryotic data. Nevertheless, the user can change the parameters. In step **4.c** a selection of Eukaryotic Metagenome-Assembled Genomes takes place. **(4.c.1)** First, the Eukaryotic Metagenome-Assembled Genomes are grouped by taxonomy using the results from **BUSCO** and **EukCC**. **(4.c.2)** Furthermore, the produced cluster are separated by species using ANI (Average Nucleotide Identity) splitter with default ANI_distance 0.95. Clusters with bootstrap > 75 are represent new species. **(4.c.3)** Finally, there is a selection of representative Metagenome-Assembled Genomes, from the clusters produced in **(4.c.2)**.

## Module 5: Relative abundace 
Module 5 is divided in 4 steps:

### 5.a: Calculation of relative abundance and genome coverage of Prokaryotic Metagenome-Assembled Genomes and construction of relative abundance and genome coverage tables
In this step the relative abundance of Prokaryotic Metagenome-Assembled Genomes is calculated and a relative abundance table is constructed (Stenetorp et al., 2012). First, **(5.a.1)** the Paired-End (PE) reads from each library are merged. Then, **(5.a.2)** the libraries are mapped to the indexed prokaryotic bins and **(5.a.3)** the hits are counted. In the next step, **(5.a.4)** a 2-column crosstable is created, with libraries as columns and prokaryotic bin representatives as rows and the percentage of the unmapped reads is calculated. After this, **(5.a.5)** the genome coverage of the mapped reads is calculated, using the BRAT results, the average number of base pairs per library and the number of base pairs and the total contigs per bin.

### 5.b: Calculation of relative abundance and genome coverage of Uncultivated Viral Genomes and construction of relative abundance and genome coverage tables 
In this step the relative abundance of Uncultivated Viral Genomes is calculated using BRAT and a relative abundance table is constructed (Stenetorp et al., 2012). First, **(5.b.1)** the Paired-End (PE) reads from each library are merged. Then, **(5.b.2)** the libraries are mapped to the indexed viral contigs and **(5.b.3)** the hits are counted. In the next step, **(5.b.4)** a 2-column crosstable is created, with libraries as columns and viral contig representatives as rows and the percentage of the unmapped reads is calculated. After this, **(5.b.5)** the genome coverage of the mapped reads is calculated, using the BRAT results, the average number of base pairs per library and the number of base pairs and the total contigs per bin.

### 5.c: Calculation of relative abundance and genome coverage of Eukaryotic Metagenome-Assembled Genomes and construction of relative abundance and genome coverage tables
In this step, the relative abundance of Eukaryotic Metagenome-Assembled Genomes is calculated a relative abundance table is constructed (Stenetorp et al., 2012). First, **(5.c.1)** the Paired-End (PE) reads from each library are merged. Then, **(5.c.2)** the libraries are mapped to the indexed eukaryotic bins and **(5.c.3)** the number of hits is counted. In the next step, **(5.c.4)** a 2-column crosstable is created, with libraries as columns and eukaryotic bin representatives as rows and the percentage of the unmapped reads is calculated. After this, **(5.c.5)** the genome coverage of the mapped reads is calculated, using the BRAT results, the average number of base pairs per library and the number of base pairs and the total contigs per bin.

### 5.d: Construction of combined relative abundance and combined genome coverage tables 
In the beginning of this step, **(5.d.1)** the combined relative abundance table of the tables constructed in **5.a**, **5.b** and **5.c** is generated. Then, **(5.d.2)** the combined genome coverage table from the genome coverage tables constructed in **5.a**, **5.b** and **5.c** is formatted.

## References
* Agarwala R, Morgulis A: BMTagger aka Best Match Tagger is for removing human reads from metagenomics datasets. In., ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/ (2010).
* Bolduc, B. Clustering Viral Genomes in iVirus. (2017) doi:10.17504/protocols.io.gwebxbe.
* Bornemann, T. L. V., Esser, S. P., Stach, T. L., Burg, T. & Probst, A. J. uBin – a manual refining tool for metagenomic bins designed for educational purposes. bioRxiv 2020.07.15.204776 (2020) doi:10.1101/2020.07.15.204776.
* Galiez, C., Siebert, M., Enault, F., Vincent, J. & Söding, J. WIsH: who is the host? Predicting prokaryotic hosts from metagenomic phage contigs. Bioinformatics 33, 3113–3114 (2017). 
* Simão, F. A., Waterhouse, R. M., Ioannidis, P., Kriventseva, E. V. & Zdobnov, E. M. BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics 31, 3210–3212 (2015).
* Holt, C. & Yandell, M. MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects. BMC Bioinformatics 12, 491 (2011).
* Alneberg, J. et al. Binning metagenomic contigs by coverage and composition. Nat Methods 11, 1144–1146 (2014).
* Jang, H. B. et al. Gene sharing networks to automate genome-based prokaryotic viral taxonomy. bioRxiv 533240 (2019) doi:10.1101/533240.
* Brown, J., Pirrung, M. & McCue, L. A. FQC Dashboard: integrates FastQC results into a web-based, interactive, and extensible FASTQ quality control tool. Bioinformatics 33, 3137–3139 (2017).
* Kang, D. D. et al. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ 7, e7359 (2019).
* Kieft, K., Zhou, Z. & Anantharaman, K. VIBRANT: automated recovery, annotation and curation of microbial viruses, and evaluation of viral community function from genomic sequences. Microbiome 8, 90 (2020).1.
* Krueger, F. FelixKrueger/TrimGalore. (2018).
* Li, D. et al. MEGAHIT v1.0: A fast and scalable metagenome assembler driven by advanced methodologies and community practices. Methods 102, 3–11 (2016).
* Nayfach, S. et al. A genomic catalog of Earth’s microbiomes. Nature Biotechnology 1–11 (2020) doi:10.1038/s41587-020-0718-6.
* Nayfach, S., Camargo, A. P., Eloe-Fadrosh, E., Roux, S. & Kyrpides, N. CheckV: assessing the quality of metagenome-assembled viral genomes. bioRxiv 2020.05.06.081778 (2020) doi:10.1101/2020.05.06.081778.
* Nurk, S., Meleshko, D., Korobeynikov, A. & Pevzner, P. MetaSPAdes: A new versatile metagenomic assembler. Genome Research 27, gr.213959.116 (2017).
* Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P. & Tyson, G. W. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Res 25, 1043–1055 (2015).1.
* Parks, D. H. et al. A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. Nat Biotechnol 36, 996–1004 (2018).
* Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P. & Parks, D. H. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. Bioinformatics 36, 1925–1927 (2020).
* Stenetorp, P. et al. BRAT: a web-based tool for NLP-assisted text annotation. in Proceedings of the Demonstrations at the 13th Conference of the European Chapter of the Association for Computational Linguistics 102–107 (Association for Computational Linguistics, 2012).
* Ren, J., Ahlgren, N. A., Lu, Y. Y., Fuhrman, J. A. & Sun, F. VirFinder: a novel k-mer based tool for identifying viral sequences from assembled metagenomic data. Microbiome 5, 69 (2017).
* Roux, S. A Viral Ecogenomics Framework To Uncover the Secrets of Nature’s “Microbe Whisperers”. mSystems 4, (2019).
* Roux, S., Enault, F., Hurwitz, B. L. & Sullivan, M. B. VirSorter: mining viral signal from microbial genomic data. PeerJ 3, e985 (2015).
* Saary, P., Mitchell, A. L. & Finn, R. D. Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC. Genome Biology 21, 244 (2020).
* Seemann, T. Prokka: rapid prokaryotic genome annotation. Bioinformatics 30, 2068–2069 (2014).
* Sieber, C. M. K. et al. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. Nat Microbiol 3, 836–843 (2018).
* Uritskiy, G. V., DiRuggiero, J. & Taylor, J. MetaWRAP-a flexible pipeline for genome-resolved metagenomic data analysis. Microbiome 6, 158 (2018).
* West, P. T., Probst, A. J., Grigoriev, I. V., Thomas, B. C. & Banfield, J. F. Genome-reconstruction for eukaryotes from complex natural microbial communities. bioRxiv 171355 (2017) doi:10.1101/171355.
* Wu, Y.-W., Simmons, B. A. & Singer, S. W. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics 32, 605–607 (2016).




 








