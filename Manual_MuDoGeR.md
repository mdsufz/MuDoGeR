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

Please note that you can click the tool's name on the text to be redirected to the tool's own documentation

## Required Metadata table

Currently, MuDoGeR v1.0 only works with paired-end sequences. Future updates will add tools to work with long-read sequencing samples.
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
## Please make sure your foward and reverse raw sequencing files end in *"_1.fastq"* and in *"_2.fastq"*, respectively.

### Naming your samples: good practices

Please pay attention to your samples' names. Avoid using spaces and special characters in your samples' names as good practices. If possible, use only letters and numbers. Special and space characters could cause malfunctioning in some of the used tools.


Following you have an usage tutorial for each module. Each **Module's final consideration** section has the description of the expected folder and files names output.

**Keep in mind that you only need one single command to run each module**

## Module 1: Pre-Processing

For running all module 1 use

```console
$ mudoger --module preprocess --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20 --metaspades -m 100
```

The available parameters for module 1 are:
* --meta metadata table as described [here](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#required-metadata-table)
* -o output directory (mandatory)
* -m  Given RAM to assembly (optional)
* -t number of threads/cores (default = 1)
* Assembler. Can be --megahit or --metaspades (default)

The -m parameter is optional, once it is calculated during step **1.b**. If you specify the amount of RAM, the result from **1.b** will be ignored. 
**Be aware that RAM is a relevant bottleneck for MAG recovery.**

### 1.a: Raw Read Quality Control  
For quality control, MuDoGeR uses the implementation present in [metaWRAP](https://github.com/bxlab/metaWRAP). Trim Galore is used for adapter and quality trimming with default PHRED quality and stringency. Future updates will allow the option to remove host contamination using different databases.


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

If you want to use the other MuDoGeR modules, you can copy the showed folder structure with your own data and use the same ```-o output/path ``` when running the other modules.

### Most relevant files from Module 1

The relevant files generated in Module 1 used for the other modules are ```assembly/final_assembly.fasta```, ```qc/final_pure_reads_1.fastq```, and ```qc/final_pure_reads_2.fastq```. Please, make sure your resulted files are specificaly named as ```final_assembly.fasta``` and ```final_pure_reads_1/2.fastq```.


## Module 2: Recovery of Prokaryotic Metagenome-Assembled Genomes (pMAGs)

This module integrates the prokariotic binning procedure implemented in [**metaWRAP**](https://github.com/bxlab/metaWRAP) using [**Metabat2**](https://peerj.com/articles/7359/), [**Maxbin2**](https://academic.oup.com/bioinformatics/article/32/4/605/1744462), and [**CONCOCT**](https://www.nature.com/articles/nmeth.3103) binners, a MuDoGeR bin dereplacation method, the [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182) taxonomic annotation, [**CheckM**](https://genome.cshlp.org/content/25/7/1043) quality estimation, [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) metagenomic gene annotation,  MuDoGeR sequence metrics and [**BBtools**](https://sourceforge.net/projects/bbmap/) sequence metrics calculation tools.

The tools that require specific databases in this module are [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182) and [**CheckM**](https://genome.cshlp.org/content/25/7/1043). Both should be ready to use after running the database-setup.sh script. See instructions [here](https://github.com/mdsufz/MuDoGeR#installation)

For running all steps in module 2 use

```console
$ mudoger --module prokaryotes --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20
```
The available parameters for module 2 are:
* --meta metadata table as described [here](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#required-metadata-table)
* -o output directory (mandatory)
* -t number of threads/cores (default = 1)

Additional modularity for this module is scheduled to happen.

**Please refer to Module 2 final consideration to understand the folder structure from Module 2.**

### 2.a: Prokaryotic sequences binning

The binnig process starts by using [**Metabat2**](https://peerj.com/articles/7359/), [**Maxbin2**](https://academic.oup.com/bioinformatics/article/32/4/605/1744462), and [**CONCOCT**](https://www.nature.com/articles/nmeth.3103) to bin the sequences from the ```final_assembly.fasta``` file.
Following, the results from all binners are used to refine bacterial bins. For bacterial bins, the refinement process uses 50% minimum completeness and 10% maximum contamination as default. For archeal bins, the refinement process uses 40% minimum completeness and 30% maximum contamination as default. The refinement process used is implemented in [**metaWRAP**](https://github.com/bxlab/metaWRAP). Finally, MuDoGeR removes redundant bins.

The most relevant files from this step are the prokaryotic bins inside the ```unique_bins/``` folder. 

### 2.b: Taxonomic classification, quality estimation, and gene annotation.

After the binning step is completed, the resulted bins are taxonomic annotated using the [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182) software and its most updated database. Following, the ```unique_bins/``` are checked for quality using the [**CheckM**](https://genome.cshlp.org/content/25/7/1043) tool. Finally, the recovered prokaryotic bins are annotated using [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517).

Inside the ```metrics/``` folder, you should have one folder for each tool. Inside ```prokka/``` you will find one folder for each bin containing the outputs from the tool.
The resulted files will be used by MuDoGeR to generate a more comprehensive report of the bins, as well as further processing. If you would like to know more about the outputs of each tool, please check their respective documentation. You can find links to the tools [here](https://github.com/mdsufz/MuDoGeR/blob/master/dependencies_description.md) or by clicking on the tool's name.

### 2.c: Sequence metrics calculation and selection of Prokaryotic MAGs.

Finally, MuDoGeR calculates some relevant metrics from the recovered bins, such as genome_size, number_of_scaffolds, largest_scaffold_size, N50, and N90. In addition, it also counts the number of annotated and unknown genes by prokka. [**BBtools**](https://sourceforge.net/projects/bbmap/) is also used to extract sequence metrics. Later, MuDoGeR merges the results from the other tools and calculates the sequence quality (completeness – 5×contamination [(Parks, 2018)](https://www.nature.com/articles/s41564-017-0012-7)). Bins with quality greater than or equal to 50 are considered MAGs and have their information summarised in the ```mags_results_summary.tsv``` file.

The ```mags_results_summary.tsv``` contains relevant annotations from the recovered MAGs. You can also check the annotated genes for each MAG by looking at the ```.tsv``` output by [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) for each MAG.

## Module 2 final considerations

By the end, MuDoGeR parses the outputs from all used tools and outputs the following comprehensive results to the ```sample_name/prokaryotes/final_outputs/``` folder. Please go to the [understand outputs](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md) to better understand MuDoGeR results.

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
	     ├── final_outputs
	     │   ├── all_bins_seq
	     │   ├── allbins_metrics_summary.tsv
	     │   ├── bins_genes_prokka_summary
	     │   ├── bins_metrics_summary
	     │   │   ├── qual_bins_checkm_summary.tsv
	     │   │   └── taxa_bins_gtdbtk_summary.tsv
	     │   ├── mags_results_summary.tsv
	     │   └── only_mags_seq
             └── metrics
                 ├── checkm_qc
		 │   ├── bins
		 │   ├── lineage.ms
		 │   ├── outputcheckm.tsv
		 │   └── storage
                 ├── genome_statistics
                 │   ├── bbtools.tsv
                 │   ├── genome_metrics.tsv
                 │   └── prok_genomes_stats.tsv
                 ├── GTDBtk_taxonomy
		 │   └── GTDB-tk_outputs
                 └── prokka
		     ├── bin.0.fa/
                     ├── bin.10.fa/
                     └── bin.11.fa/

```

## Module 3: Recovery of Uncultivated Viral Genomes (UViGs)

**Please refer to Module 3 final consideration to understand the folder structure from Module 3.**

The recovery of the UViGs module integrates the viral sequence recovery tools [**VirSorter2**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y), [**VirFinder**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0283-5), and [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0). Later, the sequences are dereplicated using [**Stampede-clustergenomes**](https://bitbucket.org/MAVERICLab/stampede-clustergenomes/src/master/). Following the putative viral contigs are analyzed and taxonomy is estimated using [**Vcontact2**](https://www.nature.com/articles/s41587-019-0100-8), and quality is estimated with [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7). MuDoGeR also uses [**WiSH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377) to estimate the viral-host pairs from the putative viral contigs. Finally, MuDoGeR compiles the outputs from the used tools and selects high-quality UViGs as defined by [Roux, S., et al.(2019)](https://www.nature.com/articles/nbt.4306) and [Nayfach, S., et al. (2021)](https://www.nature.com/articles/s41587-020-00774-7).

The tools that require specific databases in the viral module are [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0), [**WiSH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377), and [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7). All the databases should be ready to use after running the database-setup.sh script. See instructions [here](https://github.com/mdsufz/MuDoGeR#installation)

For running all stpe in module 3 use

```console
$ mudoger --module viruses --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20
```
The available parameter for module 2 are:
* --meta metadata table as described [here](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#required-metadata-table) (mandatory)
* -o output directory (mandatory)
* -t number of threads/cores (default = 1)

### 3.a: Recovery of Putative Viral Contigs

In **3.a**, the viral recovery tools  [**VirSorter2**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00990-y), [**VirFinder**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0283-5), and [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0) are applied to the assembly fasta file ```final_assembly.fasta``` created during the preprocess module. Sequences recovered with **VirFinder** with p-value > 0.01 and length < 1000 bp are removed. Later, the independent results of each tool are combined and dereplicated with [**Stampede-clustergenomes**](https://bitbucket.org/MAVERICLab/stampede-clustergenomes/src/master/) using 70% minimum coverage and 95% minimum identity.

The result files ```vibrant_filtered_data.txt```, ```virfinder_filtered_data.txt```, and ```virsorter2_filtered_data.txt``` contain the names from the sequences identified as putative viral contigs. Those are the files used during the next steps.

The final outputs from this step are inside the ```dereplication/``` folder. The file ```dereplication/uvigs_95-70.fna``` contains the viral sequences dereplicated, and the file ```dereplication/uvigs_95-70.clstr``` contains the clusters formed from all identified putative viral contigs. 


### 3.b: Taxonomic and Quality estimation of putative viral contigs

Initially, MuDoGeR uses the ```dereplication/uvigs_95-70.fna``` file in **prodigal** to produce the amino acid sequence from the putative viral contigs. These amino acid sequences are input for the [**Vcontact2**](https://www.nature.com/articles/s41587-019-0100-8) analysis and taxonomical estimation. Following, viral sequence quality is calculated using  [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7).

The main output from [**Vcontact2**](https://www.nature.com/articles/s41587-019-0100-8) used by MuDoGeR is the ```vcontact-output/genome_by_genome_overview.csv```. There you can find the estimated viral taxonomy and confidence score. For a more detailed explanation, please check the [**Vcontact2**](https://www.nature.com/articles/s41587-019-0100-8) documentation. Finally, you can find the outputs from **CheckV** inside the ```vcheck_quality``` folder. The resulted file from  [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7) is organized in ```vcheck_quality/quality_summary.tsv```.


### 3.c: Viral-Host pair estimation

The Viral-Host pair estimation process is done based on the [**WIsH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377) software. For MuDoGeR to run it automatically, please keep the defined folder structure from the Prokaryotes recovery and Viral recovery. MuDoGeR will use the results from Module 2 and Module 3 to calculate the Viral-Host potential. 

Inside the ```host_prediction/``` folder you will find the Models created by [**WIsH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377) for the prokaryotic MAGs. Inside the ```potential_host_genomes/``` you will find the MAGs recovered in Module 2, and inside the ```uvigs/``` folder, you have the fasta files for the viral contigs recovered in Module 3. Finally, you will have your main output in ```host_prediction/output_results/prediction.list```. There you will have the viral contigs name, the best provided MAGs match, the LogLikelihood, and the p-value for the match.

### 3.d: Selection of UViGs

Finally, MuDoGeR retrieves all the previously created files and creates the ```sample_name/viruses/final_outputs/``` folder. A quality filter is also applied.

The ```putative_vir_contigs_summary.tsv``` contains the compiled information from each putative viral contig (original_contig, uvig_length, provirus, proviral_length, gene_count, viral_genes, host_genes, checkv_quality, miuvig_quality, completeness, completeness_method, contamination, kmer_freq, warnings, putative_host, and likelihood). The ```uvigs_high_quality_summary.tsv``` is a subset from ```putative_vir_contigs_summary.tsv``` with only defined UViGs.

## Module 3: Final Considerations

By the end, MuDoGeR parses the outputs from all used tools and outputs the following comprehensive results to the ```sample_name/viruses/final_outputs/``` folder. Please go to the [understand outputs](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md) to better understand MuDoGeR results.

After a successful run of Module 3 you should have the following folder structure:

After step 3.a is finished, you should have the following folder structure:

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

After step 3.b is finished, you should have the following folder structure:

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

After step 3.c is finished, you should have the following folder structure:

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
     	  ├──final_outputs
	  │	   ├── only_uvigs_seq
	  │	   │    ├── sample_name_putative_viral_contig-15.fa
	  │	   │    ├── sample_name_putative_viral_contig-16.fa
	  │	   │    └── ...
	  │	   ├── putative_vir_contigs_summary.tsv
	  │	   ├── putative_vir_seq_metrics_summary
	  │	   │    ├── host_vir_pair_wish_summary.tsv
	  │	   │    ├── quality_checkv_summary.tsv
	  │	   │    ├── taxa_estimation_vir_vcontact2_summary.csv
	  │	   │    ├── vir_contigs_annotation_vibrant_summary.tsv
	  │	   │    └── vir_contigs_genbank_annotation_vibrant.tsv
	  │	   ├── putative_viral_contigs
	  │	   │    ├── sample_name_putative_viral_contig-0.fa
	  │	   │    ├── sample_name_putative_viral_contig-1.fa
	  │	   │    └── …
	  │	   ├── uvigs_high_quality_summary.tsv
	  │	   └── viral_contigs_seq_names.csv
          ├── host_prediction
          ├── investigation
          ├── taxonomy
          └── vcheck_quality

```

## Module 4: Recovery of Eukaryotic Metagenome-Assembled Bins (eMABs)

**Please refer to Module 4 final consideration to understand the folder structure from Module 4.**

The recovery of eukaryotic genomes can be trickier than the prokaryotic genome recovery. This means that we do not assume the recovered eukaryotic bins to be genomes at any moment. That is why we keep the nomenclature as eMABs instead of eMAGs. The recovery of the eMABs module integrates [**EukRep**](https://genome.cshlp.org/content/28/4/569) for selecting the eukaryotic contigs from the initial assembly, the [**CONCOCT**](https://www.nature.com/articles/nmeth.3103) binner, [**GeneMark**](https://academic.oup.com/nar/article/29/12/2607/1034721) for prediction of eukaryotic genes, [**EukCC**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02155-4) for for quality estimation from eukaryiotic sequences, [**MAKER2**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491) for gene annotation, and [**BUSCO**](https://academic.oup.com/bioinformatics/article/31/19/3210/211866) for detection of single-copy orthologous genes.

The tools that require specific databases in the eukaryotic module are [**EukCC**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02155-4), [**BUSCO**](https://academic.oup.com/bioinformatics/article/31/19/3210/211866), and [**MAKER2**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491). [**BUSCO**](https://academic.oup.com/bioinformatics/article/31/19/3210/211866) and [**EukCC**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02155-4) databases  should be ready to use after running the database-setup.sh script. See instructions [here](https://github.com/mdsufz/MuDoGeR#installation). 

Module 4 requires specific configuration that can't be done automatically. Please, make sure you follow the instructions [here](https://github.com/mdsufz/MuDoGeR/blob/master/installation/genemark_maker2_installation.md) to complete [**GeneMark**](https://academic.oup.com/nar/article/29/12/2607/1034721) and [**MAKER2**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491) installation.

For running all module 4 use:

```console
$ mudoger --module eukaryotes --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20
```

The available parameters for module 4 are:
* --meta metadata table as described [here](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#required-metadata-table) (mandatory)
* -o output directory (mandatory)
* -t number of threads/cores (default = 1)

### 4.a: Recovery and binning of Eukaryotic assemblies

The eMABs recovery starts with ```final_assembly.fasta``` as input to the [**EukRep**](https://genome.cshlp.org/content/28/4/569) tool. This splits the contigs within the assembled sequences into ```prokaryotic_contigs.fa``` and ```eukaryotic_contigs.fa```. Following, the eukaryotic contigs serve as input to the [**CONCOCT**](https://www.nature.com/articles/nmeth.3103) binner. Finally, the resulted eukaryotic bins are filtered by minimum size (by default we use 1.5 MB) and considered eMABs.

The main result of step 4.a is the eMABs fasta files located within the ```filtered_euk_bins/``` folder.

### 4.b: Completeness and contamination estimation and annotation of Eukaryotic bins

The second part of eMABs recovery starts by predicting eukaryiotic genes using [**GeneMark**](https://academic.oup.com/nar/article/29/12/2607/1034721). The input for this step are the eMABs recovered in step 4.a located within the ```filtered_euk_bins/``` folder. Following [**EukCC**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02155-4) is used to calculate quality and completeness from the eMABs recovered in 4.a.

For the annotation of eukaryotic genes, MuDoGeR feeds the [**MAKER2**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491) tool with the eukaryotic bins and the necessary files generated by [**GeneMark**](https://academic.oup.com/nar/article/29/12/2607/1034721) for each eMAB. Finally, [**BUSCO**](https://academic.oup.com/bioinformatics/article/31/19/3210/211866) receives the [**MAKER2**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491) genes as input for detecting single-copy orthologous genes (SCGs)

The results from [**GeneMark**](https://academic.oup.com/nar/article/29/12/2607/1034721) are calculated for each eMAB and you may find them in ```eukaryotes/genemarker_annotation/eMAB.fa_genemark```. The annotated genes are located in the ```genemark.gtf``` file. The model created for the eMAB, the ```gmhmm.mod``` is then used by [**MAKER2**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491).
A summary of the quality results calculated buy [**EukCC**](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02155-4) can be found in ```eukaryotes/eukcc_quality/eMAB.fa_eukcc/eukcc.csv```.
The main results from the [**MAKER2**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-491) annotation can be found for each eMAB within the folder ```eukaryotes/maker2_gene_annotation/eMAB.fa_maker2/eMAB.maker.output/```
Finally, you can find [**BUSCO**](https://academic.oup.com/bioinformatics/article/31/19/3210/211866) results in the files ```eukaryotes/euk_completeness/eMAB_busco/short_summary.specific.eukaryota_odb10.eMAB_busco.txt``` and ```eukaryotes/euk_completeness/eMAB_busco/run_eukaryota_odb10/full_table.tsv```.

## Module 4: Final Considerations

By the end, MuDoGeR parses the outputs from all used tools. Please go to the [understand outputs](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md) to better understand MuDoGeR results.

After a successful run of Module 4 you should have the following folder structure:

After step 4.a is finished, you should have the following folder structure:

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

After step 4.b is finished, you should have the following folder structure:

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

## Module 5: Relative abundance 

Module 5 of MuDoGeR maps the quality-controlled reads to the recovered pMAGs/UViGs/eMAB to calculate their abundance within the WGS samples.

Consequently, MuDoGeR uses the path to the quality-controlled reads and the path to the recovered pMAGs/UViGs/eMAB sequences.

Initially, module 5 selects the representatives' pMAGs/UViGs/eMAB based on the Average nucleotide identity (ANI) 95 distance measure using the **gOTUpick** software integrated into the MuDoGeR custom scripts. Following, MuDoGeR provides two mapping pipelines, called ```--reduced``` and ```--complete```.

When using the ```--reduced``` flag, MuDoGeR maps the representative pMAGs/UViGs/eMAB only to the WGS samples where they were recovered. When using the ```--complete``` flag, MuDoGeR maps all WGS samples to all representative pMAGs/UViGs/eMAB selected by the **gOTUpick**.  In both cases, you can use the ```--coverage``` and ``` --relative-abundance``` flags to also calculate the pMAGs/UViGs/eMABs coverage and relative abundance, respectively.

In addition, you can use the ``` --gene``` flag to map the WGS quality-controlled reads to all the [Prokka]( https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) annotated genes found in the samples' ```final_assembly.fasta``` files. You can also use the ``` --coverage```  and ``` --relative-abundance```  flags to output coverage and relative abundance values, respectively, for the annotated genes.

You can run all module 5 as follows:

```console
$ mudoger --module abundance_tables --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20 --reduced --absolute-values --coverage --relative-abundance
```

The available parameters for module 5 are:

* --meta metadata table as described [here](https://github.com/mdsufz/MuDoGeR/blob/master/Manual_MuDoGeR.md#required-metadata-table) (mandatory)
* -o output directory (mandatory)
* -t number of threads/cores (default = 1)
* --reduced (default), --complete, or --genes mapping type
* --absolute-values --coverage, and/or --relative-abundance output values

Please refer to Module 5 final consideration to understand the folder structure from Module 5.

### 5.a: Select representative pMAGs from each created OTU

The first step from module 5 is to group all recovered pMAGs within OTU groups. To do so, MuDoGeR uses the ANI 95 distance strategy implemented in the **gOTUpick** tool. Essentially, it starts by calculating the ANI distance from the pMAGs and separates them into a group according to 95% identity within the groups. Following, it uses information from the [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182) taxonomical annotation, and [**CheckM**](https://genome.cshlp.org/content/25/7/1043) quality estimation to select the highest quality pMAG from the same taxonomical class within the same ANI95 group as the group's representative MAG.
After step 5.a is finished you should find the ```mapping_results/gOTUpick_results/final_output/``` folder. Inside this folder, you will have the ```bestbins.txt``` and ```final_groups_output.csv``` files. Those files indicate the representatives' pMAGs, marked with an \*, and their groups. Those are also the files used in the next steps of the pipeline. Please go to the [understand outputs]( https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md) to better understand MuDoGeR results.

### 5.b: pMAGs/UViGs/eMABs mapping and abundance calculation

Following the MuDoGeR copies the OTUs selected in step **5.a**, the selected UViGs and eMABs recovered in module 3 and module 4, to the ```pmags_otu_mapping/otus_fasta```, ```uvigs_mapping/uvigs_fasta```, and ```euk_mabs_mapping/emabs_fasta```, respectively. 
Later, MuDoGeR uses the **Bowtie2** software to index the copied fasta files. 
If ```--coverage``` is selected, MuDoGeR calculates the pMAGs/UViGs/eMABs sizes, and the average read length of all samples to be mapped. Consequently, you should find the file ```merged_reads/avg_reads_len.tsv``` containing the average read length of all used samples. In addition, you should find the files ```genome_sizes```, ```uvigs_sizes.txt ```, and ```emabs_sizes.txt ``` inside the folder relative to the pMAGs, UViGS, and eMABs mapping, respectively. If ```--relative-abundance``` is selected, it also calculates the total reads per sample, and you should also find the file ``` merged_reads/total_reads_per_lib.tsv```.

Following this, the MuDoGeR does the actual mapping according to the ```--reduced``` or ```--complete``` strategy. You should find the folders ```map_results_ reduced``` or ```map_results_complete``` inside the folder for each domain. These folders will contain one file for each performed mapping named “mapped_bin_name”-LIB-“sample_used”.txt.
Finally, MuDoGeR calculates the absolute number of hit, coverage, and relative abundance tables and outputs the results inside the ```map_final_tables_reduced``` and ``` map_final_tables_complete```. You should have the tables named as ```map_reduced(complete)_absolute_n_hits(coverage)(relative_abundance)_table.tsv```.

Please go to the [understand outputs](https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md) to better understand MuDoGeR results.

### 5.c: Genes relative abundance calculation from the samples’ assembly

If you chose the ```--genes``` mapping type, MuDoGeR will copy and index the sample’s ```final_assembly.fasta``` file. After that, you should have the ```assembly_gene_map/raw_mapping/``` folder containing the assembly fasta file and the index outputs. Later, MuDoGeR maps the assembled sequence and you should have the ```sample_name.map.sorted.bam``` file inside the ```raw_mapping``` folder. 
Following, MuDoGeR annotates the sample’s ```final_assembly.fasta``` file using [Prokka]( https://academic.oup.com/bioinformatics/article/30/14/2068/2390517). It then converts the resulted .gff file from [Prokka]( https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) into .gtf file. You should find the results from this functional annotation step inside the ```assembly_gene_map/functional_annotation/sample_name/``` folder.
MuDoGeR then counts the mapped reads on each gene using the information present on the .gtf file and the mapped ```.map.sorted.bam``` file. The result from the counting is in the ```map_absolute_count/sample_name.count``` file. 
To calculate the gene coverage and relative abundance, MuDoGeR also calculates the average read length from the used sample and the gene length. The results from this step is found on the ```avg_reads_len.tsv``` and the ```genelength/sample_name.genelength``` respectively. 
Finally, MuDoGeR calculates the gene coverage and gene TPM normalization. The results from this step are at the ```map_coverage_norm/sample_name.cov``` and ```map_tpm_norm/sample_name.tpm``` files respectively.

Please go to the [understand outputs]( https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md) to better understand MuDoGeR results.

## Module 5: Final Considerations

By the end, MuDoGeR parses the outputs from all used tools. Please go to the [understand outputs]( https://github.com/mdsufz/MuDoGeR/blob/master/understand_main_outputs.md) to better understand MuDoGeR results.
After a successful run of Module 5 you should have the following folder structure:
After step 5.a is finished, you should have the following folder structure:
```
mapping_results
└── gOTUpick_results
   ├── ANI_OTU95
   │   └── group-1-1
   ├── ANI_distances
   │   └── group-1
   ├── all_fastani-out-1500-0.txt
   ├── allgtdb.tsv
   ├── final_output
   │   ├── bestbins.txt
   │   └── final_groups_output.csv
   ├── results
   │   ├── bestbins_metadata.tsv
   │   ├── final_groups.tsv
   │   ├── final_groups_qual.tsv
   │   ├── further_95.txt
   │   ├── groups_ANI95.tsv
   │   ├── metadata_bbtools_sorted.tsv
   │   ├── metadata_checkm_sorted.tsv
   │   ├── metadata_gtdb_taxonomy_sorted.tsv
   │   └── not_further_95.txt
   └── tax_groups
       ├── group-1
       └── unique_tax
```
After step 5.b is finished, you should have the following folder structure:
```
mapping_results
├── euk_mabs_mapping
│   ├── emabs_fasta
│   │   ├── sample_name-bin.6.1.bt2
│   │   └── …
│   ├── emabs_sizes
│   │   └── sample_name-bin.6.emab-size
│   ├── emabs_sizes.txt
│   ├── map_final_tables_complete
│   │   ├── map_complete_absolute_n_hits_list.tsv
│   │   ├── map_complete_absolute_n_hits_table.tsv
│   │   ├── map_complete_coverage_list.tsv
│   │   ├── map_complete_coverage_table.tsv
│   │   ├── map_complete_relative_abundance_list.tsv
│   │   └── map_complete_relative_abundance_table.tsv
│   ├── map_final_tables_reduced
│   │   ├── map_reduced_absolute_n_hits_list.tsv
│   │   ├── map_reduced_absolute_n_hits_table.tsv
│   │   ├── map_reduced_coverage_list.tsv
│   │   ├── map_reduced_coverage_table.tsv
│   │   ├── map_reduced_relative_abundance_list.tsv
│   │   └── map_reduced_relative_abundance_table.tsv
│   ├── map_results_complete
│   │   └── sample_name-bin.6-LIB-A.txt
│   └── map_results_reduced
│       └── sample_name-bin.6-LIB-A.txt
├── merged_reads
│   ├── sample_name.fasta
│   ├── avg_reads_len.tsv
│   └── total_reads_per_lib.tsv
├── pmags_otu_mapping
│   ├── genome_size
│   │   ├── sample_name-bin.0.genome-size
│   │   └── sample_name-bin.7.genome-size
│   ├── genomes_sizes
│   ├── map_final_tables_complete
│   │   ├── map_complete_absolute_n_hits_list.tsv
│   │   ├── map_complete_absolute_n_hits_table.tsv
│   │   ├── map_complete_coverage_list.tsv
│   │   ├── map_complete_coverage_table.tsv
│   │   ├── map_complete_relative_abundance_list.tsv
│   │   └── map_complete_relative_abundance_table.tsv
│   ├── map_final_tables_reduced
│   │   ├── map_reduced_absolute_n_hits_list.tsv
│   │   ├── map_reduced_absolute_n_hits_table.tsv
│   │   ├── map_reduced_coverage_list.tsv
│   │   ├── map_reduced_coverage_table.tsv
│   │   ├── map_reduced_relative_abundance_list.tsv
│   │   └── map_reduced_relative_abundance_table.tsv
│   ├── map_results_complete
│   │   ├── sample_name-bin.0-LIB-sample_name.txt
│   │   └── sample_name-bin.7-LIB-sample_name.txt
│   ├── map_results_reduced
│   │   ├── sample_name-bin.0-LIB-sample_name.txt
│   │   └── sample_name-bin.7-LIB-sample_name.txt
│   └── otus_fasta
│       ├── A-bin.0.1.bt2
│       ├── A-bin.7.rev.1.bt2
│       └── …
└── uvigs_mapping
    ├── map_final_tables_complete
    │   ├── map_complete_absolute_n_hits_list.tsv
    │   ├── map_complete_absolute_n_hits_table.tsv
    │   ├── map_complete_coverage_list.tsv
    │   ├── map_complete_coverage_table.tsv
    │   ├── map_complete_relative_abundance_list.tsv
    │   └── map_complete_relative_abundance_table.tsv
    ├── map_final_tables_reduced
    │   ├── map_reduced_absolute_n_hits_list.tsv
    │   ├── map_reduced_absolute_n_hits_table.tsv
    │   ├── map_reduced_coverage_list.tsv
    │   ├── map_reduced_coverage_table.tsv
    │   ├── map_reduced_relative_abundance_list.tsv
    │   └── map_reduced_relative_abundance_table.tsv
    ├── map_list_reduced
    ├── map_results_complete
    │   ├── sample_name-sample_name_putative_viral_contig-15-LIB-sample_name.txt
    │   └── …
    ├── map_results_reduced
    │   ├── sample_name-sample_name_putative_viral_contig-15-LIB-sample_name.txt
    │   └── …
    ├── uvigs_fasta
    │   ├── sample_name-sample_name_putative_viral_contig-15.1.bt2
    │   └── …
    ├── uvigs_sizes
    │   ├── sample_name-sample_name_putative_viral_contig-15.uvig-size
    │   └──…
    └── uvigs_sizes.txt
```
After step 5.c is finished, you should have the following folder structure:

```
mapping_results
└── assembly_gene_map
   ├── avg_reads_len.tsv
   ├── functional_annotation
   │   └── sample_name
   ├── genelength
   │   └── sample_name.genelength
   ├── map_absolute_count
   │   └── sample_name.count
   ├── map_coverage_norm
   │   └── sample_name.cov
   ├── map_tpm_norm
   │   └── sample_name.tpm
   └── raw_mapping
       ├── sample_name.map.sam
       └── …
```



