# DETAILED DESCRIPTION OF MuDoGeR

## Module 1: Pre-Processing 
Module 1 is divided into 3 steps:

### 1.a: Raw Read Quality Control
The **MuDoGeR** Raw Read Quality Control is based on **metaWrap** Read Quality Control (Uritskiy et al., 2018) and is responsible to prepare the raw reads for assembly and alignment. The trimming of the raw reads is based on adapted content and **PHRED** scored with the default setting of Trim-galore (Krueger, 2015). As a result, all the non-high-quality sequences are removed. The high-quality reads that remained are aligned to the genome of a potential host with bmtagger (Agarwala and Morgulis, 2010), a process that leads to the removal of host contamination and of read pairs with only one aligned read from the metagenomic data. In the final step, **FASTQC** (Brown et al., 2019) is used to assess read quality improvement.

### 1.b: Calculation of resources
In case the user is interested in a fine quality assembly, then **metaSPades** (Nurk et al., 2017) should be used. Due to the high memory requirements of **metaSPades**, prior to the assembly, it is possible to estimate the necessary resources for the assembly of the library produced in **1.a**. **(1.b.1)** This can be achieved by estimating the complexity of the forward or reversed pure reads by k-mer calculation, as they are almost identical. The k-mer is usually estimated in size 33 or 55. Both k-mer values are estimated and k-mer with sizes of 33 and 55 are combined to a single output file. **(1.b.2)** Using this file, the necessary memory that **metaSPades** needs to utilize for assembly the good-quality reads is predicted.

### 1.c: Assembly
The **MuDoGeR** assembly is based on **metaWrap** assembly and is utillized for assembling a set of metagenomic reads. Two possible tools can be chosen for the assembly: **MegaHit** (default) and **metaSPAdes**. **MegaHit** (Li et al., 2016) requires less memory and is also faster, compared to **metaSPAdes** (if the user decides to utilize **metaSPAdes** we recommend the use of **1.b**). If the user is interested in a fine quality assembly, he will use **metaSPAdes** and not **MegaHit**. Due to the higher memory requirements of metaspades prior of the assembly it is possible to estimate the resources necessary for metaSPades usage. Because of that, the usage of **metaSPAdes** reader is preferable except of the cases of large data-set.

## Module 2: Recovery of Prokaryotic Metagenome-Assembled Genomes
Module 2 is separated in 3 steps:

### 2.a: Binning of Prokaryotic Metagenome-Assembled Genomes, bin_refinement, taxonomic classification, quality estimation and annotation of Prokaryotic bins
This step is a combination of tasks included in the **metaWrap** tool (Uritskiy et al., 2018). **(2.a.1)** This step starts with the binning of the assembly data-sets with the use of three wrapped tools: **MaxBin2** (Wu et al., 2016), **metaBAT2** (Kang et al., 2019) and **CONCOCT** (Alneberg et al., 2014). This process results in the production of bin files from each method. Next, **(2.a.2)** the bin_refinement of the bins produced from each tool follows. In the bin_refinement, the minimum completeness and contamination for bacterial bins are 50% and 10% respectively, while the minimum completeness and contamination for the archaeal bins are 40% and 30% respectively. **(2.a.3)** For the taxonomic classification of the bin sets from Bacteria and Archaea produced in the bin refinement module, **GTDB-Tk** tool (Chaumeil et al., 2020) is used. Furthermore, **(2.a.4)** the quality matrix of the bins produced from bin_refinement, using **CheckM** (Parks et al., 2014). After this step, **(2.a.5)** the user can choose to filter the **CheckM** bins by quality. Bin quality is defined as completeness – 5×contamination (Parks et al., 2018). By default, the quality value is 50, but the user can choose a different value. Finally, **(2.a.6)** the prediction of functional annotation of the bins produced in the bin_refinement is achieved by **PROKKA** (Seemann, 2014), which utilizes a number of software tools.

### 2.b: Selection of Prokaryotic Metagenome-Assembled Genomes Representatives
In this step, a selection of Prokaryotic Metagenome-Assembled Genomes according Nayfach et al. (2020) takes place. It uses average nucleotide identities (ANI) distances and the threshold is 0.95 (species-level OTUs). **(2.b.1)** First, the Prokaryotic Metagenome-Assembled Genomes are grouped by taxonomy. **(2.b.2)** Furthermore, the produced cluster are separated by species using ANI (Average Nucleotide Identity) splitter with default ANI_distance 0.95. Clusters with bootstrap > 75 are represent new species. **(2.b.3)** Finally, there is a selection of representative Metagenome-Assembled Genomes, from the clusters produced in **(2.b.2)**.

### 2.c: Refinement of the selected Prokaryotic Metagenome-Assembled Genomes using U-bin
In case the user chooses to  **uBin** (Bornemann et al., 2020) refining tool can be used instead of **DAS Tool** (Sieber et al., 2018) in the bins produced in **(2.a.1)**. The users have also the option of using **uBin** to manually refined the bins produced in **2.b**. The hyperlink of the **uBin** tool can be found in the following hyperlink: ![ubin](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#ubin).

## Module 3: Recovery of Viral Metagenome-Assembled Genomes
Module 3 is divided into 2 steps: 

### **3.a**: Recovery, quality estimation, taxonomic classification and host identification of Viral Metagenome-Assembled Genomes

In the beginning of this step, three viral recovery tools (**VirFinder** (Ren et al., 2017), **VirSorter** which recovers prophages (Roux et al., 2015) and **VIBRANT** (Kieft et al., 2020)) are utilized for the identification and recovery of viral genomes from an assemble metagenomic data-set. **(3.a.1)** Each of the tools is used to recover independently the viral genomes from a given assembled data-set and saves them in separate folders. **(3.a.2)** Following that, the recovered sequences are filtered. The selection of proper **VirFinder** sequences is based on low q-value (q-value=< 0.01) and high length (length>= 1000 bp), the **VirSorter** chosen sequences are those classified into categories 1 and 2 and from **VIBRANT** the selected sequences are those of the combined assemblies from phages. **(3.a.3)** The headers of **VirSorter** and **VIBRANT** contigs are modified so they can match with those produced from the **VirFinder** results. It is important to note that in contrast with the other two tools, the **VirFinder** output recovery file contains only the headers of the assemblies. Because of that, the headers of the unique filtered sequences from each tool are extracted to a common fasta file and sorted by length. **(3.a.4)** Using the headers including in the common fasta file, the actual sequences from the assembly data-set are extracted and transferred to a new fasta file. **(3.a.5)** Next, the duplicated contigs are removed by dereplication using **Stampede-clustergenomes** tool (Bolduc and Roux, 2017). The dependencies of the dereplication step are the maximum coverage (-c) 70% and the maximum identity (-i) 95%. The dereplication is followed by **(3.a.6)**  quality estimation  of the dereplicated contigs by **CheckV** (Nayfach et al., 2020). **(3.a.7)** The taxonomic classification of the clean contigs produced by **CheckV** is done by **vContact2** (Jang et al., 2019, Bornemann et a., 2020). In case of small datasets, we can use the results of the dereplication. If the user chooses to, **(3.a.8)** the hosts of the clean contigs, or in case of small datasets the dereplicated contigs, are identified by using **WIsH tool** (Galiez et al., 2017). 

### **3.b**: Selection of Viral Metagenome-Assembled Genomes Representatives 
In **(3.b)** the user can select viral representatives. The viral representatives encompass all viruses that yielded taxonomy when using vContact2 plus those larger than 15 Kb.


## Module 4: Recovery of Eukaryotic Metagenome-Assembled Genomes
Module 4 is divided into 3 steps:

### 4.a: Recovery of Eukaryotic assemblies and production of Eukaryotic bins

In the beginning this step, **(4.a.1)** the assembled data-set is separated to Prokaryotic and Eukaryotic assemblies using **EukRep** tool (West et al., 2018). **(4.a.2)** In case the size of the Eukaryotic assembly file is >= 2.5 Mb, the Eukaryotic recovery modules can continue to the automated binning with **CONCOCT** tool. The eukaryotic bins  produced by **CONCOCT**. Finally, **(4.a.3)** the bins are filtered by size and those of size < 2.5 Mb are removed. 

### 4.b: Completeness/contamination estimation and annotation of Eukaryotic bins 

In this step, a chain of processes is followed for one of the bins produced in **4.a**: **(4.b.1)** **GeneMark-ES** is applied for gene prediction. **(4.b.2)** Also, the contamination of the bins which were kept after the filtering is estimated by using **EukCC** tool (Saary et al., 2020). **(4.b.3)** The predicted genes from **GeneMark-ES** (Saary et al., 2020) are annotated with **Maker2** (Holt et al., 2011). **(4.b.4)** The completeness of the annotated proteins is measured using **BUSCO** (Simão et al., 2015). 

### 4.c: Selection of Eukaryotic Metagenome-Assembled Genomes Representatives
Currently, a selection of Eukaryotic Metagenome-Assembled Genomes Representatives has not yet been benchmarked. Therefore, we are using standards benchmarkers with prokaryotic data. Nevertheless, the user can change the parameters. In step **4.c** a selection of Eukaryotic Metagenome-Assembled Genomes takes place. **(4.c.1)** First, the Eukaryotic Metagenome-Assembled Genomes are grouped by taxonomy using the results from **BUSCO** and **EukCC**. **(4.c.2)** Furthermore, the produced cluster are separated by species using ANI (Average Nucleotide Identity) splitter with default ANI_distance 0.95. Clusters with bootstrap > 75 are represent new species. **(4.c.3)** Finally, there is a selection of representative Metagenome-Assembled Genomes, from the clusters produced in **(4.c.2)**.

## Module 5 Relative abundace (not finished yet)
### 5.a Calculation of relative abundance of Prokaryotic Metagenome-Assembled Genomes and construction of relative abundance table 
In this step the relative abundance of Prokaryotic Metagenome-Assembled Genomes is calculated and a relative abundance table is constructed (Pontus et al., 2012, Pontus et al., 2012). First the Paired-End (PE) reads from each library are merged. Then, the libraries are mapped to the indexed prokaryotic bins and the hits are counted. In the next step, a  2-column crosstable is created, with libraries in the columns and prokaryotic bins in the rows. Finally, information of the unmapped reads is added (This part still needs to be clarified). 

### 5.b Calculation of relative abundance of Viral Metagenome-Assembled Genomes and construction of relative abundance table 
In this step the relative abundance of Viral Metagenome-Assembled Genomes is calculated using BRAT and a relative abundance table is constructed (Pontus et al., 2012, Pontus et al., 2012). First the Paired-End (PE) reads from each library are merged. Then, the libraries are mapped to the indexed viral contigs and the hits are counted. In the next step, a  2-column crosstable is created, with libraries in the columns and viral contigs in the rows. Finally, information of the unmapped reads is added (This part still needs to be clarified). 

### 5.c Calculation of relative abundance of Prokaryotic Metagenome-Assembled Genomes and construction of relative abundance table 
In this step, the relative abundance of Eukaryotic Metagenome-Assembled Genomes is calculated a relative abundance table is constructed (Pontus et al., 2012, Pontus et al., 2012). First the Paired-End (PE) reads from each library are merged. Then, the libraries are mapped to the indexed eukaryotic bins and the number of hits is counted. In the next step, a  2-column crosstable is created, with libraries in the columns and eukaryotic bins in the rows. Finally, information of the unmapped reads is added (This part still needs to be clarified). 


### 5.d Construction of relative abundance table (not done yet, no script, no information)
In this step the combined relative abundance table of the tables constructed in **5.a**, **5.b** and **5.c** is generated.

## References

* Agarwala R, Morgulis A: BMTagger aka Best Match Tagger is for removing human reads from metagenomics datasets. In., ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/, 3.101 edn: Bioconda; 2010. Accessed 15 Feb 2018.
* Benjamin Bolduc, Simon Roux 2017. Clustering Viral Genomes in iVirus. protocols.io
https://dx.doi.org/10.17504/protocols.io.gwebxbe
* Bolduc B, Jang H Bin, Doulcier G, You Z, Roux S, Sullivan MB. (2017). vConTACT: an iVirus tool to classify
double-stranded DNA viruses that infect Archaea and Bacteria. PeerJ 5: e3243.
* Bornemann, Till L.V., Sarah P. Esser, Tom L. Stach, Tim Burg, and Alexander J. Probst. “UBin – a Manual Refining Tool for Metagenomic Bins Designed for Educational Purposes.” BioRxiv, January 1, 2020, 2020.07.15.204776. https://doi.org/10.1101/2020.07.15.204776.
* Clovis Galiez, Matthias Siebert, François Enault, Jonathan Vincent, Johannes Söding, WIsH: who is the host? Predicting prokaryotic hosts from metagenomic phage contigs, Bioinformatics, Volume 33, Issue 19, 01 October 2017, Pages 3113–3114, https://doi.org/10.1093/bioinformatics/btx383.
* Ewing B, Hillier L, Wendl MC, Green P. Base-calling of automated sequencer traces using phred. I. Accuracy assessment. Genome Res. 1998 Mar;8(3):175-85. doi: 10.1101/gr.8.3.175. PMID: 9521921.
* Felipe A. Simão, Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia V. Kriventseva, Evgeny M. Zdobnov, BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs, Bioinformatics, Volume 31, Issue 19, 1 October 2015, Pages 3210–3212, https://doi.org/10.1093/bioinformatics/btv351
* Holt, C., Yandell, M. MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects. BMC Bioinformatics 12, 491 (2011). https://doi.org/10.1186/1471-2105-12-491
* Johannes Alneberg, Brynjar Smári Bjarnason, Ino de Bruijn, Melanie Schirmer, Joshua Quick, Umer Z Ijaz, Leo Lahti, Nicholas J Loman, Anders F Andersson & Christopher Quince. 2014. Binning metagenomic contigs by coverage and composition. Nature Methods, doi: 10.1038/nmeth.3103
* Jang, H. Bin, Bolduc, B., Zablocki, O., Kuhn, J., Roux, S., Adriaenssens, E., … Sullivan, M. (2019).
Gene sharing networks to automate genome-based prokaryotic viral taxonomy. BioRxiv, 533240. https://doi.org/10.1101/533240
* Joseph Brown, Meg Pirrung, Lee Ann McCue, FQC Dashboard: integrates FastQC results into a web-based, interactive, and extensible FASTQ quality control tool, Bioinformatics, Volume 33, Issue 19, 01 October 2017, Pages 3137–3139, https://doi.org/10.1093/bioinformatics/btx373 
* Kang DD, Li F, Kirton E, et al. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ. 2019;7:e7359. Published 2019 Jul 26. doi:10.7717/peerj.7359.
* Kieft, K., Zhou, Z. & Anantharaman, K. VIBRANT: automated recovery, annotation and curation of microbial viruses, and evaluation of viral community function from genomic sequences. Microbiome 8, 90 (2020). https://doi.org/10.1186/s40168-020-00867-0.
* Krueger F. Trim Galore!: a wrapper tool around Cutadapt and FastQC to
consistently apply quality and adapter trimming to FastQ files. In.,http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/, 0.4.5 edn: Bioconda; 2015. Accessed 15 Feb 2018.
* Li D, Luo R, Liu CM, Leung CM, Ting HF, Sadakane K, Yamashita H, Lam TW. MEGAHIT v1.0: A fast and scalable metagenome assembler driven by advanced methodologies and community practices. Methods. 2016 Jun 1;102:3-11. doi: 10.1016/j.ymeth.2016.02.020. Epub 2016 Mar 21. PMID: 27012178.
* Nayfach, S., Roux, S., Seshadri, R. et al. A genomic catalog of Earth’s microbiomes. Nat Biotechnol (2020). https://doi.org/10.1038/s41587-020-0718-6.
* Nayfach, Stephen, Antonio Pedro Camargo, Emiley Eloe-Fadrosh, Simon Roux, and Nikos Kyrpides. “CheckV: Assessing the Quality of Metagenome-Assembled Viral Genomes.” BioRxiv, January 1, 2020, 2020.05.06.081778. https://doi.org/10.1101/2020.05.06.081778.
* Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. metaSPAdes: a new versatile metagenomic assembler. Genome Res. 2017;27(5):824-834. doi:10.1101/gr.213959.116.
* Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2014. Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25: 1043-1055.
* Parks, D.H., Rinke, C., Chuvochina, M. et al. Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life. Nat Microbiol 2, 1533–1542 (2017). https://doi.org/10.1038/s41564-017-0012-7.
* Pierre-Alain Chaumeil, Aaron J Mussig, Philip Hugenholtz, Donovan H Parks, GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database, Bioinformatics, Volume 36, Issue 6, 15 March 2020, Pages 1925–1927, https://doi.org/10.1093/bioinformatics/btz848.
Pontus Stenetorp, Sampo Pyysalo, Goran Topić, Tomoko Ohta, Sophia Ananiadou and Jun'ichi Tsujii (2012). brat: a Web-based Tool for NLP-Assisted Text Annotation. In Proceedings of the Demonstrations Session at EACL 2012.
Pontus Stenetorp, Sampo Pyysalo, Goran Topić, Sophia Ananiadou and Akiko Aizawa. Normalisation with the Brat Rapid Annotation Tool. In Proceedings of the 5th International Symposium on Semantic Mining in Biomedicine, Zürich, Switzerland, September 2012.
* Ren, J., Ahlgren, N.A., Lu, Y.Y. et al. VirFinder: a novel k-mer based tool for identifying viral sequences from assembled metagenomic data. Microbiome 5, 69 (2017). https://doi.org/10.1186/s40168-017-0283-5.
* Roux et al. (2015), VirSorter: mining viral signal from microbial genomic data. PeerJ 3:e985;DOI 10.7717/peerj.985.
* Saary, P., Mitchell, A.L. & Finn, R.D. Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC. Genome Biol 21, 244 (2020).
https://doi.org/10.1186/s13059-020-02155-4.
* Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics. 2014 Jul 15;30(14):2068-9. doi: 10.1093/bioinformatics/btu153. Epub 2014 Mar 18. PMID: 24642063.
* Sieber CMK, Probst AJ, Sharrar A, Thomas BC, Hess M, Tringe SG, Banfield JF. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. Nat Microbiol. 2018 Jul;3(7):836-843. doi: 10.1038/s41564-018-0171-1. Epub 2018 May 28. PMID: 29807988; PMCID: PMC6786971.
* Uritskiy, G.V., DiRuggiero, J. & Taylor, J. MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis. Microbiome 6, 158 (2018). https://doi.org/10.1186/s40168-018-0541-1.
* West, Patrick T., Alexander J. Probst, Igor V. Grigoriev, Brian C. Thomas, and Jillian F. Banfield. “Genome-Reconstruction for Eukaryotes from Complex Natural Microbial Communities.” BioRxiv, January 1, 2017, 171355. https://doi.org/10.1101/171355.
* Yu-Wei Wu, Blake A. Simmons, Steven W. Singer, MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets, Bioinformatics, Volume 32, Issue 4, 15 February 2016, Pages 605–607, https://doi.org/10.1093/bioinformatics/btv638.





 








