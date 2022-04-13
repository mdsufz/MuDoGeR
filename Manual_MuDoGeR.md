# MuDoGeR 1.0 Manual

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
$ mudoger --module preprocess --meta /path/to/metadata.tsv -o /path/to/output/folder -t 20
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



## Module 3: Recovery of Uncultivated Viral Genomes
Note: Make sure that all the viral tools are correctly installed. The links for the installation can be found in the following hyperlink: ![Viral module](https://github.com/mdsufz/MuDoGeR#viral-module).

### 3.a: Recovery, quality estimation, taxonomic classification and host identification of Uncultivated Viral Genomes
In **3.a**, the viral recovery tools **VirSorter**, **VirFinder** and **VIBRANT** are applied to the assembly fasta file, for the recovery of the viral genomes contained in that. The indepentent results of each tool, combined and dereplicated with **Stampede-clustergenomes**. Also, the user can estimate the quality of the dereplicated viral contigs (**CheckV**). Also, the user can do the taxonomic classification of the the clean contigs produced by **CheckV** or in case of small data-sets the dereplicated viral contigs with **vContact2**. Furthermore, the user can choose to determine the prokaryotic host of each virus. Before running the script, it is important for the user to decide about the parameters of minimum coverege (-c) and minimum identity (-i) used in the dereplication. By default, the minimum coverage is 70% and the minimum identity 95%. However the user is free to change the dereplication parameters depending on the aims of the metagenomic analysis or the assembled dataset. 

Running **3.a**:
``` 
mudoger 3.a  -o /path/to/output/folder -f ~/path/to/assembly/file -c 70 -i 95 --g /path/to/prokaryotic_hosts/folder
```
* The `/path/to/output/folder` indicates the path to the output directory where the output folders of **3.a** will be written.
* The `/path/to/assembly/file` indicates the path to the file of the assemblies. 
* The `/path/to/prokaryotic_hosts/folder` indicates the path to the directory that contains the genomes of the possible prokaryotic hosts (optional).
* The `-c` indicates the minimum coverege.
* The `-i` indicates the minimum identity.

In the output directory five folders are present. The `initial_recovey_folder` contains the results from the independent recovery of each tool. The `dereplication_folder` includes the dereplication results using **Stampede-clustergenomes**, while in the `taxonomy_folder` and `quality_folder` the user can find the results of the taxonomic classification utlizing **vContact2** tool and quality control using **CheckV** tool, respectively. In case the user chooses to use the **WiSH** tool, a fifth folder will be present, called `wish_folder`. This folder will contain the results of the host identification analysis of either the clean contigs produced by **CheckV** or in case of small data-sets the dereplicated contigs.

```
dereplication_folder  initial_recovey_folder  quality_folder  taxonomy_folder  wish_folder 	
``` 
The `VIRAL_PARTICLES_95-70.clstr` file contains the header and the length of the contigs. For example `head -5 dereplication_folder/VIRAL_PARTICLES_95-70.clstr`:

```
>Cluster_0	NODE_6_length_9839_cov_5.151165  	   9839
>Cluster_1	NODE_7_length_8843_cov_4.647815  	   8843
>Cluster_2	NODE_9_length_8285_cov_3.135358  	   8285
>Cluster_3	NODE_16_length_6117_cov_3.296107 	   6117
>Cluster_4	NODE_20_length_5444_cov_5.065875 	   5444
```
The sequences of the contigs can be found in the `VIRAL_PARTICLES_95-70.fna` file.

The quality summary of the contigs can be found inside the `quality_folder/quality_summary.csv` file. 
Running:

```
head -5 quality_folder/quality_summary.csv 
contig_id,contig_length,genome_copies,gene_count,viral_genes,host_genes,checkv_quality,miuvig_quality,completeness,completeness_method,contamination,prophage,termini
NODE_6_length_9839_cov_5.151165,9839,1.0,16,5,0,Low-quality,Genome-fragment,16.25,AAI-based,0.0,No,No
NODE_7_length_8843_cov_4.647815,8843,1.0,10,5,0,Low-quality,Genome-fragment,15.0,AAI-based,0.0,No,No
NODE_9_length_8285_cov_3.135358,8285,1.0,14,1,0,Not-determined,Genome-fragment,NA,NA,0.0,No,No
NODE_21_length_5441_cov_6.763832,5441,1.0,8,8,0,Medium-quality,Genome-fragment,89.8,AAI-based,0.0,No,55-bp-DTR
```
In the `taxonomy_folder` the user can find the output inside the `taxonomy/v-contact_output/genome_by_genome_overview.csv` file. Running:

``` 
head -5 taxonomy/v-contact_output/genome_by_genome_overview.csv
,Genome,Order,Family,Genus,VC,VC Status,Size,VC Subcluster,VC Subcluster Size,Quality,Adj P-value,Topology Confidence Score,Genera in VC,Families in VC,Orders in VC,Genus Confidence Score
0,Achromobacter~phage~83-24,Caudovirales,Siphoviridae,Jwxvirus,0_0,Clustered,2,VC_0_0,2,0.193,0.95227493,0.1838,1,1,1,1.0
1,Achromobacter~phage~JWAlpha,Caudovirales,Podoviridae,Jwalphavirus,7_1,Clustered,11,VC_7_1,11,0.4917,1.0,0.4917,3,1,1,0.9636
2,Achromobacter~phage~JWX,Caudovirales,Siphoviridae,Jwxvirus,0_0,Clustered,2,VC_0_0,2,0.193,0.95227493,0.1838,1,1,1,1.0
3,Achromobacter~phage~phiAxp-1,Caudovirales,Siphoviridae,Unassigned,1_0,Clustered,12,VC_1_0,12,0.7098,1.0,0.7098,4,1,1,1.0
``` 

The results of the host identification are located to `wish_folder/output_results/prediction.list` file. Running:

```
head -5 wish_folder/output_results/prediction.list
"Phage"	"Best hit among provided hosts"	"LogLikelihood"	"p-value if null parameters provided"
viral-particle-260	ERR1341880_bacbin.1	-1.31937	NA
viral-particle-216	ERR1341880_bacbin.1	-1.33232	NA
viral-particle-55	LS08Hbin.1	-1.34156	NA
viral-particle-241	ERR1341880_bacbin.1	-1.29327	NA
``` 
### 3.b: Selection of Uncultivated Viral Genomes

In this step, a selection of Uncultivated Viral Genomes takes place. The Viral representatives are all viral genomes that yielded taxonomic classification with vContact2 and are larger than 15 Kb. A bash script is required for the selection.

Running **3.b**:

```
mudoger 3.b  -o /path/to/output/file -i /path/to/genome_by_genome_overview_csv -s 15
```
* The `/path/to/output/file` indicates the path to the output directory where the results of **3.b** will be written.
* The `/path/to/genome_by_genome_overview_csv` indicates the path to the input file with the results from the vcontact2 tool.
* `-s` indicates the minimum size for the filtering of viral contigs.

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

Note: Make sure that all the necessary tools are correctly installed. The links for the installation can be found in the following hyperlink: ![Relative abundance](https://github.com/mdsufz/MuDoGeR/blob/master/README.md#relative-abundance).

### 5.a Calculation of relative abundance and genome coverage of Prokaryotic Metagenome-Assembled Genomes and construction of relative abundance and genome coverage tables
For calculation of relative abundance and construction of relative abundance table for the Prokaryotic Metagenome-Assembled Genomes Representatives, the user can run:
```
mudoger 5.a -i ~/path/to/representative_bins/folder -l ~/path/to/libraries/folder -o ~/path/to/output/folder 
``` 
* The `/path/to/representative_bins/folder` indicates the path to the representative bins folder.
* The `/path/to/libraries/folder` indicates the path to the folder with the libraries of the sample. 
* The `/path/to/output/folder` indicates the path to the output folder where the results will be saved.

The output folder contains the results of the **5.a** in .sam file format. Also, inside the output folder, there is a `mappings` folder which contains the file `brat_v1.csv` which is a crosstable with libraries in the columns and bins in the rows.

Also, there is the optional step **5.a.5**, in which the user can calculate the genome coverage:

Running **5.a.5**: (it is an opinion) 

```
mudoger 5.a.5 -A ~/path/to/abundance_table -B ~/path/to/basepairs_total_contigs_per_bin -C ~/path/to/average_number_of_basepairs_per_library -o ~/path/to/output/folder 
```
* The `/path/to/abundance_table` indicates the path to the relative abundance table.
* The `/path/to/basepairs_total_contigs_per_bin` indicates the path to the table with the base pairs and the total contigs per bin. 
* The `/path/to/average_number_of_basepairs_per_library` indicates the path to the table with the average number of base pairs per library.
* The `/path/to/output/folder` indicates the path to the output folder where the results will be saved.

The result is the `brats_abs_cov.csv` file.

### 5.b Calculation of relative abundance and genome coverage of Uncultivated Viral Genomes and construction of relative abundance and genome coverage tables
For calculation of relative abundance and construction of relative abundance table for the Viral Metagenome-Assembled Genomes Representatives, the user can run:
```
mudoger 5.b -i ~/path/to/representative_contigs/folder -l ~/path/to/libraries/folder -o ~/path/to/output/folder 
``` 
* The `/path/to/representative_contigs/folder` indicates the path to the representative contigs folder.
* The `/path/to/libraries/folder` indicates the path to the folder with the libraries of the sample. 
* The `/path/to/output/folder` indicates the path to the output folder where the results will be saved.

The output folder contains the results of the **5.b** in .sam file format. Also, inside the output folder, there is a `mappings` folder which contains the file `brat_v1.csv` which is a crosstable with libraries in the columns and contigs in the rows.

Also, there is the optional step **5.b.5**, in which the user can calculate the genome coverage:

Running **5.b.5**: (it is an opinion)

```
mudoger 5.b.5 -A ~/path/to/abundance_table -B ~/path/to/basepairs_total_contigs_per_bin -C ~/path/to/average_number_of_basepairs_per_library -o ~/path/to/output/folder 
```
* The `/path/to/abundance_table` indicates the path to the relative abundance table.
* The `/path/to/basepairs_total_contigs_per_bin` indicates the path to the table with the base pairs and the total contigs per bin. 
* The `/path/to/average_number_of_basepairs_per_library` indicates the path to the table with the average number of base pairs per library.
* The `/path/to/output/folder` indicates the path to the output folder where the results will be saved.

The result is the `brats_abs_cov.csv` file.

### 5.c Calculation of relative abundance and genome coverage of Eukaryotic Metagenome-Assembled Genomes and construction of relative abundance and genome coverage  tables
For calculation of relative abundance and construction of relative abundance table for the Eukaryotic Metagenome-Assembled Genomes Representatives, the user can run:
```
mudoger 5.c -i ~/path/to/representative_bins/folder -l ~/path/to/libraries/folder -o ~/path/to/output/folder 
``` 
* The `/path/to/representative_bins/folder` indicates the path to the representative bins folder.
* The `/path/to/libraries/folder` indicates the path to the folder with the libraries of the sample. 
* The `/path/to/output/folder` indicates the path to the output folder where the results will be saved.

The output folder contains the results of the **5.c** in .sam file format. Also, inside the output folder, there is a `mappings` folder which contains the file `brat_v1.csv` which is a crosstable with libraries in the columns and bins in the rows.

Also, there is the optional step **5.c.5**, in which the user can calculate the genome coverage:

Running **5.c.5**:

```
mudoger 5.c.5 -A ~/path/to/abundance_table -B ~/path/to/basepairs_total_contigs_per_bin -C ~/path/to/average_number_of_basepairs_per_library -o ~/path/to/output/folder 
```
* The `/path/to/abundance_table` indicates the path to the relative abundance table.
* The `/path/to/basepairs_total_contigs_per_bin` indicates the path to the table with the base pairs and the total contigs per bin. 
* The `/path/to/average_number_of_basepairs_per_library` indicates the path to the table with the average number of base pairs per library.
* The `/path/to/output/folder` indicates the path to the output folder where the results will be saved.

The result is the `brats_abs_cov.csv` file. 

### 5.d Construction of combined relative abundance and combined genome coverage tables
(to be tested).

```
mudoger 5.d.1  -o ~/path/to/output/folder -A ~/path/to/prok_relative_abundance_table -B ~/path/to/viral_relative_abundance_table -C ~/path/to/euk_relative_abundance_table  
``` 

* The `/path/to/prok_relative_abundance_table` indicates the path to the prokaryotic relative abundance table.
* The `/path/to/viral_relative_abundance_table` indicates the path to the viral relative abundance table.
* The `/path/to/euk_relative_abundance_table` indicates the path to the eukaryotic relative abundance table.
* The `/path/to/output/folder` indicates the path to the output folder where the results will be saved.

The final table is found to the output folder as: `final_brat.csv` 

Also, there is the optional step 5.d.2, in which the user can create the combined genome coverage table:

Running **5.d.2**:

```
mudoger 5.d.2 -o ~/path/to/output/folder -A ~/path/to/prok_genome_covarage_table -B ~/path/to/viral_genome_covarage_table -C ~/path/to/euk_genome_covarage_table  
``` 

* The `/path/to/prok_genome_covarage_table` indicates the path to the prokaryotic coverage abundance table.
* The `/path/to/viral_genome_covarage_table` indicates the path to the viral coverage abundance table.
* The `/path/to/euk_genome_covarage_table` indicates the path to the eukaryotic coverage abundance table.
* The `/path/to/output/folder` indicates the path to the output folder where the results will be saved.

The final table is found to the output folder as:`final_brats_abs_cov.csv` file.


## 2.b: Selection of Prokaryotic Metagenome-Assembled Genomes Representatives
In this step, the user can pick representative Prokaryotic Metagenome-Assembled Genomes. The threshold for ANI clustering is set by default at 95 but the user has the option to change this value.

Running 2.b:
```
mudoger 2.b -i ~/path/to/gtdb_taxonomy/file -b ~/path/to/bbtools/file -m ~/path/to/bins(mags)/folder -o ~/path/to/output/folder -t 95
``` 
* The `/path/to/gtdb_taxonomy/file` indicates the path to the taxonomy file, generated by **GTDB-Tk**. The user should choose its .tsv form.
* The `/path/to/bbtools/file` indicates the path to the file with bbtools.
* The `/path/to/bins(mags)/folder`  indicates the path to the the prokaryotic bins (mags) folder.
* The `/path/to/output/folder` indicates the path to the output folder where the resulted files will be saved.
* The `-t 95` indicates the threshold for ANI clustering.

Inside the output folder the user can find the `bestbins.csv` file that contains the unique taxonomic bins and the `bins_to_brats.txt` file with the selected bins that were chosens as representatives and will be used for the Bin Relative Abundance Table (BRAT) calculation.


