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

MuDoGeR version 1.0 was designed to be an easy-to-use genome recovery tool. Therefore, we created a setup procedure that requires little user input. Consequently, one important aspect of MuDoGeR usage is the output folder architecture. The folder and files names and architecture are important for the pipeline to work smoothly. If you want to prepare your data using other tools, and later use MuDoGeR, please keep the folder structure and file naming according to the MuDoGeR requirements. 

Following you have a usage tutorial for each module. Each **Module's final consideration** section has the description of the expected folder and files names output.

## Module 1: Pre-Processing

For running module 1 use

```console
$ mudoger --module preprocess --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20 --metaspades -m 100
```
The available parameter for module 1 are:
* -o output directory
* -m  Given RAM to assembly
* -t number of threads/cores
* Assembler. Can be --megahit [default] or --metaspades

--meta metadata table containing the sample names and its raw reads path in tsv format as follows:
```console
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

### 1.a: Raw Read Quality Control  
For quality control, MuDoGeR uses the implementation present in [metaWRAP](https://github.com/bxlab/metaWRAP). The quality control procedure currently applied is to trim raw reads based on adapted content and PHRED scored with the default setting of Trim-galore. Future updates will allow the option to remove host contamination using different databases.


The output directory of the read quality control module (***sample_name/qc***) contains:
```
final_pure_reads_1.fastq    pre-QC_report
final_pure_reads_2.fastq    post-QC_report 
```

The `final_pure_reads` files contain the sequences of the trimmed and decontaminated reads. The `pre-QC_report` and `post-QC_report` folders include the html reports for the reads before and after the read quality control. 

Raw reads:

![](https://github.com/mdsufz/MuDoGeR/blob/master/Read_QC_before_trimming.png)


Reads after QC:

![](https://github.com/mdsufz/MuDoGeR/blob/master/Read_QC_after_trimming.png)


### 1.b: Calculation of resources 
After quality control, MuDoGeR estimates the required RAM for assembly.

Initially, it counts the 33 and 55 unique k-mers from the quality-controlled reads (final_pure_reads_1.fastq) in the pre-processed reads (forward or reversed).
 
The k-mer count is used by a linear regression model to predict the amount of memory necessary for assembling the reads with **metaSPAdes**. 
Inside the calculation of resources output folder (***sample_name/khmer***) the user can find the following files:
```
final_prediction.tsv  input_for_predictR.tsv  kmer-33  kmer-55  metaspades_prediction.tsv
```
The ```final_prediction.tsv``` file contains the final RAM estimation in MB for assembling the sample with **metaspades**


### 1.c: Assembly 
There are two possible tools for assembling data: **MegaHiT** and **metaSPAdes**. Both tools are considered reliable. **MegaHiT** uses lower memory and is faster compared to **metaSPAdes**, but **metaSPAdes** tends to produced longer contigs. In case of very large data sets, the user may use **MegaHiT**.

The assembly output folder (***sample_name/assembly***) should contain the following:

```
assembly_report.html  final_assembly.fasta  megahit  metaspades  QUAST_out
```

The ```final_assembly.fasta``` is the assembled sequences that are going to be used in the other modules.

![](https://github.com/mdsufz/MuDoGeR/blob/master/Assembly.png)

## Module 1 final considerations

Please, be aware of the folder structure expected as a result of module 1.

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

If you want to use the other MuDoGeR modules, you can copy the showed folder structure with your own data and use the same *** -o output/path *** when running the other modules. 
The relevant files generated in Module 1 used for the other modules are ***assembly/final_assembly.fasta***, ***qc/final_pure_reads_1.fastq***, and ***qc/final_pure_reads_2.fastq***. Please, make sure your resulted files are specificaly as final_assembly.fasta and final_pure_reads_1/2.fastq.


## Module 2: Recovery of Prokaryotic Metagenome-Assembled Genomes (MAGs)

This module integrates the prokariotic binning procedure implemented in **metaWrap** using **Metabat2**, **Maxbin2**, and **CONCOCT** binners, a MuDoGeR bin dereplacation method, the **GTDB-tk** taxonomic annotation, **CheckM** quality estimation, **Prokka** metagenomic gene annotation,  MuDoGeR sequence metrics and **BBtools** sequence metrics calculation tools.

The tools that require specific databases in this module are **GTDB-tk** and **CheckM**. Both should be ready to use after running the database-setup.sh script. See instructions [here](https://github.com/JotaKas/MuDoGeR#installation)

For running module 2 use

```console
$ mudoger --module prokaryotes --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20
```
Additional modularity for this module is scheduled to happen.

## 2.a: Prokaryotic sequences binning

The binnig process starts by using **Metabat2**, **Maxbin2**, and **CONCOCT** to bin the sequences from the ***final_assembly.fasta*** file.
Following, the results from all binners are used to refine bacterial bins. For bacterial bins, the refinement process uses 50% minimum completeness and 10% maximum contamination as default. For archeal bins, the refinement process uses 40% minimum completeness and 30% maximum contamination as default. The refinement process used is implemented in **metaWrap**. Finally, MuDoGeR removes redundant bins. 

After a successful run of prokaryotic sequence binning you should have the following folder structure:
```
sample_name
    └── prokaryotes
        └── binning
            ├── initial-binning
            ├── refinement-arc
            ├── refinement-bac
            └── unique_bins

```
Inside the ***unique_bins*** folder you should have all unique prokaryotic bins found in your sample.


## 2.b: Taxonomic classification, quality estimation, and gene annotation.

After the binning step is completed, the resulted bins are taxonomic annotated using the **GTDB-tk** software and its most updated database. Following, the unique_bins are checked for quality using the **CheckM** tool. Finally, the recovered prokaryotic bins are annotated using **Prokka**.

After successfully running step 2.b, you should have the following folder structure:

```
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
             └── metrics
                 ├── checkm_qc
                 │   ├── bins
                 │   ├── lineage.ms
                 │   ├── outputcheckm.tsv
                 │   └── storage
                 ├── GTDBtk_taxonomy
                 │   ├── align
                 │   ├── classify
                 │   ├── gtdbtk.ar122.markers_summary.tsv
                 │   ├── gtdbtk.bac120.classify.tree 
                 │   ├── gtdbtk.bac120.filtered.tsv 
                 │   ├── gtdbtk.bac120.markers_summary.tsv 
                 │   ├── gtdbtk.bac120.msa.fasta 
                 │   ├── gtdbtk.bac120.summary.tsv 
                 │   ├── gtdbtk.bac120.user_msa.fasta 
                 │   ├── gtdbtk.failed_genomes.tsv 
                 │   ├── gtdbtk.log
                 │   ├── gtdbtk.translation_table_summary.tsv 
                 │   ├── gtdbtk.warnings.log
                 │   └── identify
                 └── prokka
                     ├── bin.0.fa/
                     ├── bin.10.fa/
                     └── bin.11.fa/

```

Inside the ***metrics*** folder, you should have one folder for each tool. Inside ***prokka*** you will find one folder for each bin containing the outputs from the tool.
The resulted files will be used by MuDoGeR to generate a more comprehensive report of the bins, as well as further processing. If you would like to know more about the outputs of each tool, please check their respective documentation.


## 2.c: Sequence metrics calculation and selection of Prokaryotic MAGs.

Finally, MuDoGeR calculates some relevant metrics from the recovered bins, such as genome_size, number_of_scaffolds, largest_scaffold_size, N50, and N90. In addition, it also counts the number of annotated and unknown genes by prokka. **BBtools** is also used to extract sequence metrics. Later, MuDoGeR merges the results from the other tools and calculates the sequence quality (completeness – 5×contamination (Parks, 2018)). Bins with quality greater than or equal to 50 are considered MAGs and have their information summarised in the ***MAGS_results.tsv*** file.

After successfully running step 2.c, you should have the following folder structure:

```
sample_name
     └── prokaryotes
             ├── binning
             │   ├── initial-binning
             │   ├── refinement-arc
             │   ├── refinement-bac
             │   └── unique_bins
             ├── MAGS_results.tsv
             └── metrics
                 ├── checkm_qc
                 ├── genome_statistics
                 │   ├── bbtools.tsv
                 │   ├── genome_metrics.tsv
                 │   └── prok_genomes_stats.tsv
                 ├── GTDBtk_taxonomy
                 └── prokka

```
The ***MAGS_results.tsv*** contains relevant annotations from the recovered MAGs. You can also check the annotated genes for each MAG by looking at the ***.gtf*** output by ***Prokka*** for each MAG.



## Module 3: Recovery of Uncultivated Viral Genomes (Uvigs)

The recovery of the Uvigs module integrates the viral sequence recovery tools **VirSorter**, **VirFinder**, and **VIBRANT**. Later, the sequences are dereplicated using **Stampede-clustergenomes**. Following the Uvigs are analyzed and taxonomically estimated using ***Vcontact2**, and quality is estimated with **CheckV**. MuDoGeR also uses **WiSH** to predict the prokaryotic hosts from the Uvigs. Finally, MuDoGeR compiles the outputs from the used tools and selects high-quality Uvigs.

The tools that require specific databases in the viral module are **VIBRANT**, **WiSH**, and **CheckV**. All the databases should be ready to use after running the database-setup.sh script. See instructions [here](https://github.com/JotaKas/MuDoGeR#installation)

For running module 3 use

```console
$ mudoger --module viruses --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20
```

### 3.a: Recovery of Uncultivated Viral Genomes

In **3.a**, the viral recovery tools **VirSorter**, **VirFinder** and **VIBRANT** are applied to the assembly fasta file ***final_assembly.fasta*** created during the preprocess module. Sequences recovered with **VirFinder** with p-value > 0.01 and length < 1000 bp are removed. Later, the independent results of each tool are combined and dereplicated with **Stampede-clustergenomes** using 70% minimum coverage and 95% minimum identity.


After successfully running step 3.a, you should have the following folder structure:

```
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

The files ***vibrant_filtered_data.txt***, ***virfinder_filtered_data.txt***, and ***virsorter2_filtered_data.txt*** contains the names from the sequences identified as Uvigs. Those are the files used during the next steps.

The final outputs from this step are inside the ***dereplication*** folder. The file ***dereplication/uvigs_95-70.fna*** contains the viral sequences dereplicated, and the file file ***dereplication/uvigs_95-70.clstr*** contains the clusters formed from all identified Uvigs. 


### 3.b: Taxonomic and Quality estimation of Uvigs

Initially, MuDoGeR uses the ***dereplication/uvigs_95-70.fna*** file in **prodigal** to produce the amino acid sequence from the Uvigs. These amino acid sequences are input for the **vContact** analysis and taxonomical estimation. Following, viral sequence quality is calculated using **CheckV**.

After successfully running step 3.b, you should have the following folder structure:

```
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

The main output from **vContact** used by MuDoGeR is the ***vcontact-output/genome_by_genome_overview.csv***. There you can find the estimated viral taxonomy and confidence score. For a more detailed explanation, please check the **vContact** documentation. Finally, you can find the outputs from **CheckV** inside the ***vcheck_quality*** folder. The resulted file from **CheckV** is organized in ***vcheck_quality/quality_summary.tsv***



### 3.c: Host identification of the dereplicated Uvigs

The Viral-Host identification process is done based on the **WiSH** software. For MuDoGeR to run it automatically, please keep the defined folder structure from the Prokaryotes recovery and Viral recovery. MuDoGeR will use the results from Module 2 and Module 3 to calculate the Viral-Host potential. 

After successfully running step 3.c, you should have the following folder structure:

```
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
Inside the ***host_prediction*** folder you will find the Models created by **WiSH** for the prokaryotic MAGs. Inside the ***potential_host_genomes*** you will find the MAGs recovered in Module 2, and inside the ***uvigs*** folder, you have the fasta files for the viral sequences recovered in Module 3. Finally, you will have your main output in ***host_prediction/output_results/prediction.list***. There you will have the Uvigs name, the best provided MAGs match, the LogLikelihood, and the p-value for the match.


### 3.d: Selection of Uvigs

Finally, MuDoGeR retrieves all the previously created files and creates a summary table for the recovered Uvigs. A quality filter is also applied.

After successfully running step 3.d, you should have the following folder structure:

```
sample_name
     └── viruses
          ├── host_prediction
          ├── investigation
          ├── taxonomy
          ├── Uvigs_high_quality.tsv
          ├── vcheck_quality
          └── viruses_summary.tsv

```

The ***viruses_summary.tsv*** contains the compiled information from each Uvig (original_contig, uvig_length, provirus, proviral_length, gene_count, viral_genes, host_genes, checkv_quality, miuvig_quality, completeness, completeness_method, contamination, kmer_freq, warnings, putative_host, and likelihood). The ***Uvigs_high_quality.tsv*** is a subset from ***viruses_summary.tsv*** with only high quality Uvigs.



## Module 4: Recovery of Eukaryotic Metagenome-Assembled Genomes 

Note: Make sure that all the eukaryotic tools are correctly installed. The links for the installation can be found in the following hyperlink: ![Eukaryotic module](https://github.com/mdsufz/MuDoGeR#eukaryotic-module).

### 4.a: Recovery of Eukaryotic assemblies and production of Eukaryotic bins 
In **4.a**, the **EukRep** separates the eukaryotic from the prokaryotic assemblies and then eukaryotic bins are produced by **CONCOCT**. The bins are filtered by size. Bins with size < 2.5 Mb are removed.

Running **4.a**:

```
mudoger 4.a -f ~/path/to/assembly/file --prokarya /path/to/prokaryotic/folder -o /path/to/output/folder -1 ~/path/to/final_pure_reads_1.fastq -2 ~/path/to/final_pure_reads_2.fastq 
```
* The `/path/to/assembly/file` indicates the path to the file of the assemblies. 
* The `/path/to/output/folder` indicates the path to the output directory where the output folders of the **4.a** will be saved.
* The `/path/to/prokaryotic/folder` indicates the path to a directory where the prokaryotic assemblies will be saved after the separation of eukaryotic and prokaryotic assemblies with **EukRep**. 
* The `/path/to/final_pure_reads_1.fastq` indicates the path to the file of the forward clean reads. 
* The `/path/to/final_pure_reads_2.fastq` indicates the path to the file of the reversed clean reads. 

In the output of the first step the user can find `euk_concoct_bins` folder which contains the eukaryotic bins after the filtering.

### 4.b: Completeness/contamination estimation and annotation of Eukaryotic bins 
In **4.b**, the completeness and contamination of the Eukaryotic bins produced in **4.a** are estimated. Additionally, the annotation of these bins is taking place. **4.b** starts with **GeneMark-ES** tool, used for the gene prediction in the Eukaryotic bins. As input, the user can use any of the bins produced in **4.a**. The rest of the tools used in **4.b** are **EukCC** (bin contamination estimation), **MAKER2** (annotation of the bin) and **BUSCO** (bin completeness estimation). 

Running **4.b**:

```
mudoger 4.b -f ~/path/to/concoct/bin/fasta/file -o ~/path/to/output/folder 
```
* The `/path/to/concoct/bin/fasta/file` indicates the path to the bin file.
* The `/path/to/output/folder` indicates the path to the output directory where the output folders of **4.b** will be written.

After the end of the second step the output folder contains the results from **MAKER2**, **BUSCO** and **EukCC** tools:

The results of the **EukCC** tool are located in the `eukcc/eukcc.csv` file. The `eukcc` directory is found inside the initial output directory.

The results of the **MAKER2** tool are located in the `maker/euk-ebin.maker.output/OUTPUT.all.maker.genemark.transcripts.fasta` file which contains the names and the sequences of the annotated proteins of the predicted genes . The maker directory is found  inside the `genemark` directory. 

The results of the **BUSCO** tool are located in the file `maker/busco/full_table_fbusco.csv`. The busco directory is found inside the `maker` directory.

### 4.c: Selection of Eukaryotic Metagenome-Assembled Genomes Representatives (test)
In this step, the user can pick representative Eukaryotic Metagenome-Assembled Genomes. The threshold for ANI clustering is set by default at 95 but the user has the option to change this value.

Running 4.b:
```
mudoger 4.c -i ~/path/to/busco_results/file -b ~/path/to/eukcc_results/file -m ~/path/to/bins(mags)/folder -o ~/path/to/output/folder -t 95
``` 
* The `/path/to/busco_results/file` indicates the path to the file with the results from **BUSCO**. The user should use its .tsv form
* The `/path/to/eukcc_results/file` indicates the path to the results of **EukCC**.
* The `/path/to/bins(mags)/folder`  indicates the path to the the eukaryotic bins (mags) folder.
* The `/path/to/output/folder` indicates  the path to the output folder where the resulted files will be saved.
* The `-t 95` indicates the threshold for ANI clustering.

Inside the output folder the user can find the `bestbins.csv` that contains the unique taxonomic bins and the `bins_to_brats.txt` file with the bins for the Bin Relative Abundance Table (BRAT) calculation.


## Module 5: Relative abundance 

UNDER DEVELOPMENT

