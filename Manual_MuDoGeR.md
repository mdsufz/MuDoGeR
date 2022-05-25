# MuDoGeR Manual v1.0

```console
				███    ███ ██    ██ ██████   ██████   ██████  ███████ ██████  
				████  ████ ██    ██ ██   ██ ██    ██ ██       ██      ██   ██ 
				██ ████ ██ ██    ██ ██   ██ ██    ██ ██   ███ █████   ██████  
				██  ██  ██ ██    ██ ██   ██ ██    ██ ██    ██ ██      ██   ██ 
				██      ██  ██████  ██████   ██████   ██████  ███████ ██   ██ 
						Multi-Domain Genome Recovery
							Version 1.0.0
```

MuDoGeR version 1.0 was designed to be an easy-to-use genome recovery tool. Therefore, we created a setup procedure that requires little user input. Consequently, one **important** aspect of MuDoGeR usage is the **output folder architecture**. The file names and folder architecture are important for the pipeline to work smoothly. If you want to prepare your data using other tools, and later use MuDoGeR, please keep the folder structure and file naming according to the MuDoGeR requirements. Please, check the outputs and folders created during the successful execution of each module in the following sections.

## Required Metadata table

Currently, MuDoGeR v1.0 only works with pair-end sequences. Future updates will add tools to work with long-read sequencing samples.
The pipeline requires, as input, a metadata table in tsv (tab-separated values) format containing the samples to be processed and the path in your computer to its raw sequence reads. The metadata file should have the sample name and the path to the forward reads file from the sample in one line, followed by the same sample name and the path to the reverse reads from the sample. An example metadata file is as follows:

```console
#Show the content of the metadata.tsv file
$ cat metadata.tsv

EA_ERX4593008   /path/to/EA_ERX4593008/raw_reads_1.fastq
EA_ERX4593008   /path/to/EA_ERX4593008/raw_reads_2.fastq
EA_ERX4593009   /path/to/EA_ERX4593009/raw_reads_1.fastq
EA_ERX4593009   /path/to/EA_ERX4593009/raw_reads_2.fastq
EA_ERX4593010   /path/to/EA_ERX4593010/raw_reads_1.fastq
EA_ERX4593010   /path/to/EA_ERX4593010/raw_reads_2.fastq
EA_ERX4593011   /path/to/EA_ERX4593011/raw_reads_1.fastq
EA_ERX4593011   /path/to/EA_ERX4593011/raw_reads_2.fastq

```


Following you have an usage tutorial for each module. Each **Module's final consideration** section has the description of the expected folder and files names output.

**Keep in mind that you only need one single command to run each module**

## Module 1: Pre-Processing

For running all module 1 use

```console
$ mudoger --module preprocess --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20 --metaspades -m 100
```

The available parameter for module 1 are:
* --meta metadata table as described [here](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#required-metadata-table)
* -o output directory (mandatory)
* -m  Given RAM to assembly (optional)
* -t number of threads/cores (default = 1)
* Assembler. Can be --megahit or --metaspades (default)

The -m parameter is optional, once it is calculated during step **1.b**. If you specify the amount of RAM, the result from **1.b** will be ignored. 
Be aware that RAM is a relevant bottleneck for MAG recovery.

### 1.a: Raw Read Quality Control  
For quality control, MuDoGeR uses the implementation present in [metaWRAP](https://github.com/bxlab/metaWRAP). The quality control procedure currently applied is to trim raw reads based on adapted content and PHRED scored with the default setting of Trim-galore. Future updates will allow the option to remove host contamination using different databases.


The output directory of the read quality control module (```sample_name/qc```) contains:
```console
final_pure_reads_1.fastq    pre-QC_report
final_pure_reads_2.fastq    post-QC_report 
```

The `final_pure_reads_` files contain the sequences of the trimmed and decontaminated reads. The `pre-QC_report` and `post-QC_report` folders include the html reports for the reads before and after the read quality control. 


### 1.b: Calculation of resources 
After quality control, MuDoGeR estimates the required RAM for assembly.

Initially, it counts the 33 and 55 unique k-mers from the quality-controlled reads (final_pure_reads_1.fastq) in the pre-processed reads (forward or reversed).
 
The k-mer count is used by a linear regression model to predict the amount of memory necessary for assembling the reads with [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824). 
Inside the calculation of resources output folder (```sample_name/khmer```) the user can find the following files:
```console
final_prediction.tsv  input_for_predictR.tsv  kmer-33  kmer-55  metaspades_prediction.tsv
```
The ```final_prediction.tsv``` file contains the final RAM estimation in MB for assembling the sample with **metaspades**


### 1.c: Assembly 
There are two possible tools for assembling data: [**MegaHiT**](https://academic.oup.com/bioinformatics/article/31/10/1674/177884) and [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824). Both tools are considered reliable. [**MegaHiT**](https://academic.oup.com/bioinformatics/article/31/10/1674/177884) uses lower memory and is faster compared to [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824), but [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824) tends to produced longer contigs. In case of very large data sets, the user may use **MegaHiT**.

The assembly output folder (```sample_name/assembly```) should contain the following:

```console
assembly_report.html  final_assembly.fasta  megahit  metaspades  QUAST_out
```

The ```final_assembly.fasta``` is the assembled sequences that are going to be used in the other modules.

![](https://github.com/mdsufz/MuDoGeR/blob/master/Assembly.png)

## Module 1 final considerations

Please, be aware of the folder structure expected as a result of module 1.
After a successful run of module 1 you should have the following folder structure:

```console
output/path/
         └── sample_name
             ├── assembly
             │   ├── assembly_report.html
             │   ├── final_assembly.fasta
             │   ├── megahit
             │   ├── metaspades
             │   └── QUAST_out
             ├── khmer
             │   ├── final_prediction.tsv
             │   ├── input_for_predictR.tsv
             │   ├── kmer-33
             │   ├── kmer-55
             │   └── metaspades_prediction.tsv
             └── qc
                 ├── final_pure_reads_1.fastq
                 ├── final_pure_reads_2.fastq
                 ├── post-QC_report
                 └── pre-QC_report

```

If you want to use the other MuDoGeR modules, you can copy the showed folder structure with your own data and use the same ``` -o output/path ``` when running the other modules.

### Most relevant files from Module 1

The relevant files generated in Module 1 used for the other modules are ```assembly/final_assembly.fasta```, ```qc/final_pure_reads_1.fastq```, and ```qc/final_pure_reads_2.fastq```. Please, make sure your resulted files are specificaly named as ```final_assembly.fasta``` and ```final_pure_reads_1/2.fastq```.


## Module 2: Recovery of Prokaryotic Metagenome-Assembled Genomes (MAGs)

This module integrates the prokariotic binning procedure implemented in [**metaWRAP**](https://github.com/bxlab/metaWRAP) using [**Metabat2**](https://peerj.com/articles/7359/), [**Maxbin2**](https://academic.oup.com/bioinformatics/article/32/4/605/1744462), and [**CONCOCT**](https://www.nature.com/articles/nmeth.3103) binners, a MuDoGeR bin dereplacation method, the [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182) taxonomic annotation, [**CheckM**](https://genome.cshlp.org/content/25/7/1043) quality estimation, [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) metagenomic gene annotation,  MuDoGeR sequence metrics and [**BBtools**](https://sourceforge.net/projects/bbmap/) sequence metrics calculation tools.

The tools that require specific databases in this module are [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182) and [**CheckM**](https://genome.cshlp.org/content/25/7/1043). Both should be ready to use after running the database-setup.sh script. See instructions [here](https://github.com/mdsufz/MuDoGeR#installation)

For running all steps in module 2 use

```console
$ mudoger --module prokaryotes --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20
```
The available parameter for module 2 are:
* --meta metadata table as described [here](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#required-metadata-table)
* -o output directory (mandatory)
* -t number of threads/cores (default = 1)

Additional modularity for this module is scheduled to happen.

### 2.a: Prokaryotic sequences binning

The binnig process starts by using [**Metabat2**](https://peerj.com/articles/7359/), [**Maxbin2**](https://academic.oup.com/bioinformatics/article/32/4/605/1744462), and [**CONCOCT**](https://www.nature.com/articles/nmeth.3103) to bin the sequences from the ```final_assembly.fasta``` file.
Following, the results from all binners are used to refine bacterial bins. For bacterial bins, the refinement process uses 50% minimum completeness and 10% maximum contamination as default. For archeal bins, the refinement process uses 40% minimum completeness and 30% maximum contamination as default. The refinement process used is implemented in [**metaWRAP**](https://github.com/bxlab/metaWRAP). Finally, MuDoGeR removes redundant bins.

The most relevant files from this step are the prokaryotic bins inside the ```unique_bins/``` folder


### 2.b: Taxonomic classification, quality estimation, and gene annotation.

After the binning step is completed, the resulted bins are taxonomic annotated using the [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182) software and its most updated database. Following, the ```unique_bins/``` are checked for quality using the [**CheckM**](https://genome.cshlp.org/content/25/7/1043) tool. Finally, the recovered prokaryotic bins are annotated using [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517).

Inside the ```metrics/``` folder, you should have one folder for each tool. Inside ```prokka/``` you will find one folder for each bin containing the outputs from the tool.
The resulted files will be used by MuDoGeR to generate a more comprehensive report of the bins, as well as further processing. If you would like to know more about the outputs of each tool, please check their respective documentation. You can find links to the tools [here](https://github.com/mdsufz/MuDoGeR/blob/master/dependencies_description.md).


### 2.c: Sequence metrics calculation and selection of Prokaryotic MAGs.

Finally, MuDoGeR calculates some relevant metrics from the recovered bins, such as genome_size, number_of_scaffolds, largest_scaffold_size, N50, and N90. In addition, it also counts the number of annotated and unknown genes by prokka. [**BBtools**](https://sourceforge.net/projects/bbmap/) is also used to extract sequence metrics. Later, MuDoGeR merges the results from the other tools and calculates the sequence quality (completeness – 5×contamination [(Parks, 2018)](https://www.nature.com/articles/s41564-017-0012-7)). Bins with quality greater than or equal to 50 are considered MAGs and have their information summarised in the ```MAGS_results.tsv``` file.

The ```MAGS_results.tsv``` contains relevant annotations from the recovered MAGs. You can also check the annotated genes for each MAG by looking at the ```.gtf``` output by ```PROKKA``` for each MAG.

## Module 2 final considerations

After a successful run of Module 2 you should have the following folder structure:

```console
sample_name
     └── prokaryotes
             ├── binning
             │   ├── initial-binning
             │   ├── refinement-arc
             │   ├── refinement-bac
	     │   └── unique_bins
             │       ├── bin.0.fa
             │       ├── bin.10.fa
             │       └── bin.11.fa
             ├── MAGS_results.tsv
             └── metrics
                 ├── checkm_qc
                 ├── genome_statistics
                 │   ├── bbtools.tsv
                 │   ├── genome_metrics.tsv
                 │   └── prok_genomes_stats.tsv
                 ├── GTDBtk_taxonomy
                 └── prokka
		     ├── bin.0.fa/
                     ├── bin.10.fa/
                     └── bin.11.fa/

```

## Module 3: Recovery of Uncultivated Viral Genomes (UViGs)

The recovery of the UViGs module integrates the viral sequence recovery tools [**VirSorter**](https://peerj.com/articles/985/), [**VirFinder**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0283-5), and [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0). Later, the sequences are dereplicated using [**Stampede-clustergenomes**](https://bitbucket.org/MAVERICLab/stampede-clustergenomes/src/master/). Following the potential viral contigs are analyzed and taxonomically estimated using [**Vcontact2**](https://www.nature.com/articles/s41587-019-0100-8), and quality is estimated with [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7). MuDoGeR also uses [**WiSH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377) to estimate the viral-host pairs from the potential viral contigs. Finally, MuDoGeR compiles the outputs from the used tools and selects high-quality UViGs as defined by [Roux, S., et al.(2019)](https://www.nature.com/articles/nbt.4306) and [Nayfach, S., et al. (2021)](https://www.nature.com/articles/s41587-020-00774-7).

The tools that require specific databases in the viral module are [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0), [**WiSH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377), and [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7). All the databases should be ready to use after running the database-setup.sh script. See instructions [here](https://github.com/mdsufz/MuDoGeR#installation)

For running all module 3 use

```console
$ mudoger --module viruses --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20
```
The available parameter for module 2 are:
* --meta metadata table as described [here](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#required-metadata-table) (mandatory)
* -o output directory (mandatory)
* -t number of threads/cores (default = 1)
* 
### 3.a: Recovery of Potential Viral Contigs

In **3.a**, the viral recovery tools  [**VirSorter**](https://peerj.com/articles/985/), [**VirFinder**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0283-5), and [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0) are applied to the assembly fasta file ```final_assembly.fasta``` created during the preprocess module. Sequences recovered with **VirFinder** with p-value > 0.01 and length < 1000 bp are removed. Later, the independent results of each tool are combined and dereplicated with [**Stampede-clustergenomes**](https://bitbucket.org/MAVERICLab/stampede-clustergenomes/src/master/) using 70% minimum coverage and 95% minimum identity.

The files result files ```vibrant_filtered_data.txt```, ```virfinder_filtered_data.txt```, and ```virsorter2_filtered_data.txt``` contains the names from the sequences identified as potential viral contigs. Those are the files used during the next steps.

The final outputs from this step are inside the ```dereplication/``` folder. The file ```dereplication/uvigs_95-70.fna``` contains the viral sequences dereplicated, and the file file ```dereplication/uvigs_95-70.clstr``` contains the clusters formed from all identified potential viral contigs. 


### 3.b: Taxonomic and Quality estimation of potential viral contigs

Initially, MuDoGeR uses the ```dereplication/uvigs_95-70.fna``` file in **prodigal** to produce the amino acid sequence from the potential viral contigs. These amino acid sequences are input for the [**Vcontact2**](https://www.nature.com/articles/s41587-019-0100-8) analysis and taxonomical estimation. Following, viral sequence quality is calculated using  [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7).

The main output from [**Vcontact2**](https://www.nature.com/articles/s41587-019-0100-8) used by MuDoGeR is the ```vcontact-output/genome_by_genome_overview.csv```. There you can find the estimated viral taxonomy and confidence score. For a more detailed explanation, please check the [**Vcontact2**](https://www.nature.com/articles/s41587-019-0100-8) documentation. Finally, you can find the outputs from **CheckV** inside the ```vcheck_quality``` folder. The resulted file from  [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7) is organized in ```vcheck_quality/quality_summary.tsv```.


### 3.c: Host identification of the dereplicated Potential Viral contigs

The Viral-Host pair estimation process is done based on the [**WiSH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377) software. For MuDoGeR to run it automatically, please keep the defined folder structure from the Prokaryotes recovery and Viral recovery. MuDoGeR will use the results from Module 2 and Module 3 to calculate the Viral-Host potential. 

Inside the ```host_prediction/``` folder you will find the Models created by **WiSH** for the prokaryotic MAGs. Inside the ```potential_host_genomes/``` you will find the MAGs recovered in Module 2, and inside the ```uvigs/``` folder, you have the fasta files for the viral contigs recovered in Module 3. Finally, you will have your main output in ```host_prediction/output_results/prediction.list```. There you will have the viral contigs name, the best provided MAGs match, the LogLikelihood, and the p-value for the match.

### 3.d: Selection of UViGs

Finally, MuDoGeR retrieves all the previously created files and creates a summary table for the recovered viral contigss. A quality filter is also applied.

The ```viruses_summary.tsv``` contains the compiled information from each potential viral contig (original_contig, uvig_length, provirus, proviral_length, gene_count, viral_genes, host_genes, checkv_quality, miuvig_quality, completeness, completeness_method, contamination, kmer_freq, warnings, putative_host, and likelihood). The ```Uvigs_high_quality.tsv``` is a subset from ```viruses_summary.tsv``` with only defined UViGs.

## Module 3: Final Considerations

After successfully running step 3.a, you should have the following folder structure:

```console
sample_name
     └── viruses
          └── investigation
             ├── dereplication
             │   ├── uvigs_95-70.clstr
             │   ├── uvigs_95-70.fna
             │   ├── uvigs-cover.csv
             │   ├── uvigs.fa
             │   ├── uvigs_mapping.txt
             │   ├── uvigs-nucmer.out.coords
             │   ├── uvigs-nucmer.out.delta
             │   └── viral_unique_contigs
             ├── vibrant
             │   └── VIBRANT_final_assembly
             ├── vibrant_filtered_data.txt
             ├── virfinder
             │   └── virfinder_output.tsv
             ├── virfinder_filtered_data.txt
             ├── virsorter
             │   ├── config.yaml
             │   ├── final-viral-boundary.tsv
             │   ├── final-viral-combined.fa
             │   ├── final-viral-score.tsv
             │   ├── iter-0
             │   └── log
             └── virsorter2_filtered_data.txt

```

After successfully running step 3.b, you should have the following folder structure:

```console
sample_name
     └── viruses
           ├── investigation
           │   └── CONTENT_FROM_THE_Uvigs_RECOVERY_SHOWN_BEFORE
           ├── taxonomy
           │   ├── AUX-1
           │   ├── AUX-2
           │   ├── vcontact-output
           │   │   ├── c1.clusters
           │   │   ├── c1.ntw
           │   │   ├── genome_by_genome_overview.csv
           │   │   ├── merged_df.csv
           │   │   ├── merged.dmnd
           │   │   ├── merged.faa
           │   │   ├── merged.self-diamond.tab
           │   │   ├── merged.self-diamond.tab.abc
           │   │   ├── merged.self-diamond.tab.mci
           │   │   ├── merged.self-diamond.tab_mcl20.clusters
           │   │   ├── merged.self-diamond.tab_mcxload.tab
           │   │   ├── modules_mcl_5.0.clusters
           │   │   ├── modules_mcl_5.0_modules.pandas
           │   │   ├── modules_mcl_5.0_pcs.pandas
           │   │   ├── modules.ntwk
           │   │   ├── sig1.0_mcl2.0_clusters.csv
           │   │   ├── sig1.0_mcl2.0_contigs.csv
           │   │   ├── sig1.0_mcl2.0_modsig1.0_modmcl5.0_minshared3_link_mod_cluster.csv
           │   │   ├── sig1.0_mcl5.0_minshared3_modules.csv
           │   │   ├── vConTACT_contigs.csv
           │   │   ├── vConTACT_pcs.csv
           │   │   ├── vConTACT_profiles.csv
           │   │   ├── vConTACT_proteins.csv
           │   │   └── viral_cluster_overview.csv
           │   ├── viral_genomes.faa
           │   ├── viral_genomes_g2g.csv
           │   └── viral_genomes.genes
           └── vcheck_quality
              ├── complete_genomes.tsv
              ├── completeness.tsv
              ├── contamination.tsv
              ├── proviruses.fna
              ├── quality_summary.tsv
              ├── tmp
              │   └──checkv_only_files 
              └── viruses.fna

```

After successfully running step 3.c, you should have the following folder structure:

```console
sample_name
     └── viruses
           ├── host_prediction
           │   ├── modelDir
           │   │   ├── sample_name-bin.0.mm
           │   │   └── sample_name-bin.1.mm
           │   ├── nullmodels
           │   │   ├── computeNullParameters.R
           │   │   ├── llikelihood.matrix
           │   │   └── nullParameters.tsv
           │   ├── output_results
           │   │   ├── llikelihood.matrix
           │   │   └── prediction.list
           │   ├── potential_host_genomes
           │   │   ├── sample_name-bin.0.fa
           │   │   └── sample_name-bin.1.fa
           │   └── uvigs
           │       ├── sample_name_uvig-0.fa
           │       └── sample_name_uvig-1000.fa
           ├── investigation
           │   └── CONTENT_FROM_THE_Uvigs_RECOVERY_SHOWN_BEFORE
           ├── taxonomy
           │   └── CONTENT_FROM_THE_Uvigs_TAXONOMY_SHOWN_BEFORE
           └── vcheck_quality
                └── CONTENT_FROM_THE_Uvigs_CheckV_SHOWN_BEFORE

```

After successfully running step 3.d, you should have the following folder structure:

```console
sample_name
     └── viruses
          ├── host_prediction
          ├── investigation
          ├── taxonomy
          ├── Uvigs_high_quality.tsv
          ├── vcheck_quality
          └── viruses_summary.tsv

```

## Module 4: Recovery of Eukaryotic Metagenome-Assembled Genomes (eMABs)

The recovery of the eMABs module integrates **EukRep** for selecting the eukaryotic contigs from the initial assembly, the **CONCOCT** binner, **GeneMark** for prediction of eukaryotic genes, **EukCC** for for quality estimation from eukaryiotic sequences, **MAKER2** for gene annotation, and **BUSCO** for detection of single-copy orthologous genes.

The tools that require specific databases in the eukaryotic module are **EukCC**, **BUSCO**, and **MAKER2**. **BUSCO** and **EukCC** databases  should be ready to use after running the database-setup.sh script. See instructions [here](https://github.com/mdsufz/MuDoGeR#installation). 

Module 4 requires specific configuration that can't be done automatically. Please, make sure you follow the instructions [here](https://github.com/mdsufz/MuDoGeR/blob/master/installation/genemark_maker2_installation.md) to complete **GeneMark** and **MAKER2** installation.

For running module 4 use:

```console
$ mudoger --module eukaryotes --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20
```

### 4.a: Recovery and binning of Eukaryotic assemblies

The eMABs recovery starts with ```final_assembly.fasta``` as input to the **EukRep** tool. This splits the contigs within the assembled sequences into ```prokaryotic_contigs.fa``` and ```eukaryotic_contigs.fa```. Following, the eukaryotic contigs serve as input to the **CONCOCT** binner. Finally, the resulted eukaryotic bins are filtered by minimum size (by default we use 1.5 MB) and considered eMABs.

After successfully running step 4.a, you should have the following folder structure:

```console
sample_name
     └── eukaryotes
		├── eukaryotes_bins
		│   ├── bin.0.fa
		│   ├── bin.10.fa
		│   ├── bin.11.fa
		│   └── bin.8.fa
		├── eukaryotic_contigs.fa
		├── filtered_euk_bins
		│   └── bin.8.fa
		└── prokaryotic_contigs.fa

```

The main result of step 4.a.1 is the eMABs fasta files located within the ```filtered_euk_bins/``` folder.

### 4.b: Completeness and contamination estimation and annotation of Eukaryotic bins

The second part of eMABs recovery starts by predicting eukaryiotic genes using **GeneMark**. The input for this step are the eMABs recovered in step 4.a located within the ```filtered_euk_bins/``` folder. Following **EukCC** is used to calculate quality and completeness from the eMABs recovered in 4.a.

For the annotation of eukaryotic genes, MuDoGeR feeds the **MAKER2** tool with the eukaryotic bins and the necessary files generated by **GeneMark** for each eMAB. Finally, **BUSCO** receives the **MAKER2** genes as input for detecting single-copy orthologous genes (SCGs)

After successfully running step 4.b, you should have the following folder structure:

```console
sample_name
     └── eukaryotes
		├── eukaryotes_bins
		│   ├── bin.0.fa
		│   ├── bin.10.fa
		│   ├── bin.11.fa
		│   └── bin.8.fa
		├── eukaryotic_contigs.fa
		├── eukcc_quality
		│   └── bin.8.fa_eukcc
		│       ├── eukcc.csv
		│       ├── eukcc.log
		│       └── savestate.json.gz
		├── euk_completeness
		│   ├── bin.8_busco
		│   │   ├── logs
		│   │   ├── run_eukaryota_odb10
		│   │   ├── short_summary.specific.eukaryota_odb10.bin.8_busco.json
		│   │   └── short_summary.specific.eukaryota_odb10.bin.8_busco.txt
		│   └── busco_downloads
		│       └── file_versions.tsv
		├── filtered_euk_bins
		│   └── bin.8.fa
		├── genemarker_annotation
		│   └── bin.8.fa_genemark
		│       ├── bin.8.fa
		│       ├── data
		│       ├── genemark.gtf
		│       ├── gmes.log
		│       ├── gmhmm.mod
		│       ├── info
		│       ├── output
		│       ├── run
		│       └── run.cfg
		├── maker2_gene_annotation
		│   └── bin.8.fa_maker2
		│       ├── bin.8.fa
		│       ├── bin.8.maker.output
		│       ├── maker_bopts.ctl
		│       ├── maker_evm.ctl
		│       ├── maker_exe.ctl
		│       └── maker_opts.ctl
		└── prokaryotic_contigs.fa


```

The results from **GeneMark** are calculated for each eMAB and you may find them in ```eukaryotes/genemarker_annotation/eMAB.fa_genemark```. The annotated genes are located in the ```genemark.gtf``` file. The model created for the eMAB, the ```gmhmm.mod``` is then used by **MAKER2**.
A summary of the quality results calculated buy **EukCC** can be found in ```eukaryotes/eukcc_quality/eMAB.fa_eukcc/eukcc.csv```.
The main results from the **MAKER2** annotation can be found for each eMAB within the folder ```eukaryotes/maker2_gene_annotation/eMAB.fa_maker2/eMAB.maker.output/```
Finally, you can find **BUSCO** results in the files ```eukaryotes/euk_completeness/eMAB_busco/short_summary.specific.eukaryota_odb10.eMAB_busco.txt``` and ```eukaryotes/euk_completeness/eMAB_busco/run_eukaryota_odb10/full_table.tsv```.


## Module 5: Relative abundance 

UNDER DEVELOPMENT

Module 5 of MuDoGeR maps the quality-controlled reads to the recovered MAGs/UViGs/eMAB to calculate their abundance within the WGS samples.

Consequently, MuDoGeR requires as input the path to the quality-controlled reads and the path to the recovered MAGs and UViGs sequences.

MuDoGeR provides three mapping pipelines, called --reduced , --complete , and --genes

When using the --complete flag, MuDoGeR maps the WGS reads from a sample to the recovered MAGs from the same WGS sample. Therefore, the mapping occurs only within the sample. This flag outputs an abundance table for each MAG and UViGs found in the sample.

When using the --reduced flag, MuDoGeR collects the MAGs recovered from all WGS samples provided in the --meta metadata table and groups them based on ANI 95 distance measure. Later, MuDoGeR selects the highest quality MAG within the group and maps it to the WGS reads from the samples where the MAGs from the group were recovered. This flag outputs an abundance table with samples in columns and representative MAGs and UViGs as rows. 

In addition, you can use the --gene flag to map the WGS quality-controlled reads to all the prokka annotated genes found in all MAGs in that sample.

MuDoGeR has 3 output formats for calculating abundance: --absolute-values, --coverage, and --relative-abundance. They can be used simultaneously.

You can run module 5 as follows:

mudoger --module abundance_tables --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20 --reduced --absolute-values --coverage --relative-abundance

After a succesul run of module 5, you should have the following folder structure:
```console
output/folder
	└── mapping_results
	    ├── all_bins
	    ├── all_metrics
	    └── gOTUpick_results
	    └── abundance_results
```
Inside the ```gOTUpick_results``` folder, you should find the list of MAGs selected for the mapping if the ```--reduced``` flag is used.
Inside the ```abundance_results``` folder, you should have the list of MAGs abundance per sample if you selected the ```--complete``` flag or the abundance tables calculated if you selected the ```--reduced``` flag. The output files are named according to the selected output calculation method (--absolute-values, --coverage, and --relative-abundance). If you select the --gene flag, you should find the gene abundance tables for each MAG within the prokaryotes folder for each sample. The gene abundance tables should be under a new folder named ```gene_abundaces``` inside the ```prokaryotes``` folder.


