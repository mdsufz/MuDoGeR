# aniSplitter.R

### Overview

aniSplitter calculates clusters based on a user defined Average Nucleotide Identity (ANI) threshold.

### Installation

Download aniSplitter.R script from Github
```sh
$ git clone https://github.com/felipeborim789/aniSplitter
```


### Requirements

R (https://cran.r-project.org/mirrors.html)

R Packages:
```R
install.packages("optparse")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("dendextend")
install.packages("cluster")
install.packages("reshape2")
install.packages("fpc")
```

### Usage

```sh
$ ./aniSplitter.R -d [WORK DIRECTORY] -f [FASTANI OUTPUT FILE NAME] -t [TAXONOMY FILE NAME] -a [ANI THRESHOLD]
```

### Example

```sh
$ ./aniSplitter.R -d /home/foo/bar -f fastani-output.txt -t bin_taxonomy.txt -a 95
```

### Citation

aniSplitter.R uses two functions (read.ANI and ANI.dendrogram) adapted from https://github.com/lmc297/bactaxR, thus please also cite:

Carroll, Laura M., Martin Wiedmann, Jasna Kovac. 2020. "Proposal of a Taxonomic Nomenclature for the Bacillus cereus Group Which Reconciles Genomic Definitions of Bacterial Species with Clinical and Industrial Phenotypes." mBio 11(1): e00034-20; DOI: 10.1128/mBio.00034-20
