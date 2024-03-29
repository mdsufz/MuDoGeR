## Dependencies
### Metawrap

* **Metawrap** uses python version of 2.7.14.

Tool link: https://github.com/bxlab/metaWRAP

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

The Databases can be installed by following the instructons of the following hyperlink: [Databases_metaWrap](https://github.com/bxlab/metaWRAP/blob/master/installation/database_installation.md).

### GTDB-Tk

Tool link: https://academic.oup.com/bioinformatics/article/36/6/1925/5626182

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

Tool link: https://genome.cshlp.org/content/25/7/1043

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

### PROKKA

Tool link: https://academic.oup.com/bioinformatics/article/30/14/2068/2390517

### BBTools

Tool link: https://sourceforge.net/projects/bbmap/

### VirSorter2

Tool link: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y

| Dependencies | 
|---|
| Software HMMER of version 3.1b2 | 
| Software MCL |  
| Software Metagene Annotator |
| Software MUSCLE | 
| Software BLAST+ | 
| Software DIAMOND |
| Package gcc/4/8.1-3 | 

* The software used in the **VirSorter** can be found also in the following hyperlink: [Software_VirSorter](https://github.com/simroux/VirSorter#dependencies).

### CheckV

Tool link: https://www.nature.com/articles/s41587-020-00774-7

* For **CheckV** to work the installation of Anaconda3/5.3.0 with python version of 3.6 is required.

| Dependencies | Version  |  
|---|---|
| Software DIAMOND | 2.0.4 | 
| Software HMMER | 3.3 | 
| Software Prodigal | 2.6.3 |  

* More information about **CheckV** dependencies can be found by following the webpage: https://bitbucket.org/berkeleylab/checkv/src/master/.

### vContact2

Tool link: https://www.nature.com/articles/s41587-019-0100-8

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

Tool link: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0

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

* More information about **VIBRANT** dependencies can be found in the following hyperlink: [VIBRANT_Dependencies](https://github.com/AnantharamanLab/VIBRANT#requirements-).

### Stampede-clustergenomes

Tool link: https://bitbucket.org/MAVERICLab/stampede-clustergenomes/src/master/

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

Tool link: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0283-5

* **VirFinder** uses R/3.4.4. Also the user need to have OpenMPI/2.1.2 library installed.    

| Dependencies | 
|---|
| R package glmnet |  
| R package Rcpp |  
| R package qvalue |   
| Package GCC/6.4.0-2.28 | 
| Package OpenMPI/2.1.2 | 

### WIsH

Tool link: https://academic.oup.com/bioinformatics/article/33/19/3113/3964377

* **WIsH** requires OpenMP library installed.

### EukRep

Tool link: https://genome.cshlp.org/content/28/4/569

* **EukRep** uses python of version 3.6.4.
 
| Dependencies | 
|---|
| GCC/6.4.0-2.28 |

### MAKER2

Tool link: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491

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

* As SNAP requires some more dependencies, they can be found in the following hyperlink: [Compiling_snap](https://github.com/KorfLab/SNAP#compiling).
* As RepeatModeler(1.0.4)  requires some more dependencies,they can be found in the following hyperlink: [RepeatModeler_dependencies](https://github.com/Dfam-consortium/TETools).
* As RepeatMasker (4.0.5) requires some more dependencies, they can be found in the following link: http://www.repeatmasker.org/RepeatMasker/.
* As Augustus requires some more dependencies, they can be found in the following hyperlink: [Dependencies_Augustus](https://github.com/Gaius-Augustus/Augustus#install-dependencies).
* As BEDtools requires some more dependencies, they can be found in the following link: https://bedtools.readthedocs.io/en/latest/content/installation.html.
* About Busco, the user can check the following hyperlink:  [BUSCO_Manual](https://github.com/WenchaoLin/BUSCO-Mod).
* More information about **MAKER2** dependencies can be found in the following hyperlink: [Maker2_Software prerequisites](https://github.com/wuying1984/MAKER2_PM_genome_annotation#software-prerequisites).

### BUSCO

Tool link: https://academic.oup.com/bioinformatics/article/31/19/3210/211866

* **BUSCO** uses python of version 3.7.4.

| Dependencies |  
|---|
| Package GCC/6.4.0-2.28 | 
| Package OpenMPI/2.1.2 |  
| foss/2019b |  
| NCBI BLAST+ |
| HMMER 3.1b2 |
| Augustus 3.0.x |

* As Augustus requires some more dependencies, they can be found in the following hyperlink: [Dependencies_Augustus](https://github.com/Gaius-Augustus/Augustus#install-dependencies).
* More information about **BUSCO** dependencies can be found in the following hyperlink:[BUSCO_Dependencies](https://github.com/WenchaoLin/BUSCO-Mod#setup).

### EukCC

Tool link: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02155-4

* The **EukCC** uses python of version 3.7.4.

| Dependencies for **EukCC** |  
|---|
| Package GCCcore/8.3.0 | 
| Package foss/2019b  |  
| Perl/5.30.0 |
| Anaconda2/5.3.0 |
| Package foss-2019b |

* More information about **EukCC** dependencies can be found by following the webpage: https://eukcc.readthedocs.io/en/latest/install.html#install-via-conda.

### GeneMark-ES

Tool link: https://academic.oup.com/nar/article/29/12/2607/1034721
...

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


### Bowtie

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



