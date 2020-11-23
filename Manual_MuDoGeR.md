
# Module 1: Pre-Processing module
## 1.1: Raw read QC module 
Note: For the removal of human contamination, the user will need the bmtagger hg38 index to remove the human or use another host genome for the filtering  against with the `-x` option as it can be found in the **MetaWrap** installation instructions. 

For running the raw read QC module:

``` 
mudoger read_qc module --skip-bmtagger -1 /path/to/raw_reads_1.fasta -2 /path/to/raw_reads_2.fasta  -t 24 -o /path/to/pure_reads/output/directory
```
* The `/path/to/raw_reads_1.fasta` indicates the path to file of the forward reads of the used library.
* The `/path/to/raw_reads_2.fasta` indicates the path to file of the reversed reads of the used library.
* The `/path/to/pure_reads/output/directory` indicates the path to the output directory where the clean reads of the library will be saved.
* The `-t` indicates the number of threads to be used for read quality control.

The output directory of the read quality control module contains:
```
final_pure_reads_1.fastq    pre-QC_report
final_pure_reads_2.fastq    post-QC_report 
```

The `final_pure_reads` files contain the sequences of the trimmed and decontaminated reads. The `pre-QC_report` and `post-QC_report` folders include the html reports for the reads before and after the read quality control. 


Raw reads:

![](https://github.com/mdsufz/MuDoGeR/blob/master/Read_QC_before_trimming.png)


Reads after read QC:

![](https://github.com/mdsufz/MuDoGeR/blob/master/Read_QC_after_trimming.png)



## 1.2: Resources calculation
Before the assembly, it is possible to calculate unique k-mers in the pre-processed reads(forward or reverse). The size of the k-mer that has to be investigated is usually 33 or 55. Both values have to be calculated. As the result of this task are the same using both forward and reverse reads, the user does not have to re-do the task for both of them. 

For resource calculation the user can run:

``` 
Resources module -i /path/to/final_pure_reads_1.fastq -l /path/to/output/folder 
```
 * The `/path/to/final_pure_reads_1.fasta` indicates the path to the file of the forward clean reads. 
 * The `/path/to/output/folder` indicates the path to the folder with resource calculation results.
 
The k-mer results is used as function by a linear regression model which will give the amount of memory necessary for assembly of reads by **MetaSpades**. 
Inside the output folder the user can find the  `metaspades_prediction.csv` file which has the amount of memory that **MetaSpades** utilizes for the assembly of those reads.


## 1.3: Assembly module
The reads are assembled with utilization of **MetaSpades** option flag:
```
metawrap assembly -1 /path/to/final_pure_reads_1.fasta -2 path/to//final_pure_reads_2.fastq -m 200 -t 96 --use-metaspades -o /path/to/assembled_reads/output/directory 
```

* The `/path/to/final_pure_reads_1.fasta` indicates the path to the file of the forward clean reads. 
* The `/path/to/final_pure_reads_2.fastq` indicates the path to the file of the reverse clean reads.
* The `/path/to/assembled_reads/output/directory` indicates the path to the output directory where the assemblies will be saved.
* The `-m` the amount of memory in gigabytes that the assembly process needs. 
* The `-t` indicates the number of threads to be used for this process.

After the end of the assembly process, inside the output directory the user can find the folder `Assembly_folder`. Inside this folder is the assembly file called `final_assembly.fasta` and the `assembly_report.html` file with the QUAST assembly report of the assembly module. 

![](https://github.com/mdsufz/MuDoGeR/blob/master/Assembly.png)

Using `grep > Assembly_output/assembly.fasta | head -5`, the user can see the top five longer sequences.

``` 
>NODE_1_length_369741_cov_7.638564
>NODE_2_length_360965_cov_1208.183270
>NODE_3_length_278391_cov_5.902381
>NODE_4_length_257954_cov_1138.617195
>NODE_5_length_248688_cov_1110.129452
```

# Module 2: Recovery of Prokaryotic Metagenome-Assembled Genomes
Note: Make sure that all the databases and programms required for the MetaWrap run are downloaded.The links for the installation can be found in the installation module of the ![README](https://github.com/mdsufz/MuDoGeR/blob/master/README.md) file.

## 2.1 
The run of the prokaryotic module leads to the recovery of prokaryotic genomes from the assembly dataset by utillizing the MetaWrap tool. The script of the prokaryotic module combines the algorithms of every MetaWrap module. The important parameters of minimum completion (-c) and maximum contamination(-x) for the CheckM quaity control are settled by default to 50% and 10% respectively for the recovery of bacterial contigs, while for archaeal contigs recovery are settled by default to 40% and 30% respectively.

Run the prokayotic module with using MetaWrap:
``` 
mudoger prokaryotic module -o /path/to/metawrap/output/directory -f ~/path/to/assembly/file -1 ~/path/to/final_pure_reads_1.fastq -2 ~/path/to/final_pure_reads_2.fastq --checkM_filter -q 50
```

* The `/path/to/metawrap/output/directory` indicates the path to the output directory where the output folders of the prokaryotic module will be written.
* The `/path/to/assembly/file` indicates the path to the file of the assemblies. 
* The `/path/to/final_pure_reads_1.fastq` indicates the path to the file of the forward clean reads. 
* The `path/to/final_pure_reads_2.fastq` indicates the path to the file of the reverse clean reads. 
* The `-q` indicates lower limit of quality for the filtering for the optional step of quality control (optional).

In the final output folder the user can find:

```
Annotation_Archaea   checkm_archaea   archaea_output_tax_dir    
Annotation_Bacteria  checkm_bacteria  bacteria_output_tax_dir           
```

Inside each tax_dir directory the user can find the classification of the refined bins. As an example, running `cat tax_out_dir_bacteria/classify/intermediate_results/gtdbtk.bac120.classification_pplacer.tsv| head -5`:

``` 
bin.3	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Shewanellaceae;g__Shewanella;s__
bin.5	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__;s__
bin.4	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Rhodocyclaceae;g__Azoarcus_C;s__
bin.1	d__Bacteria;p__Desulfobacterota;c__Desulfuromonadia;o__Geobacterales;f__Geobacteraceae;g__Geobacter;s__
bin.2	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Dysgonomonadaceae;g__Dysgonomonas;s__
``` 

Inside each of the checkm directories there is the `bin` directory that contains their fasta. Also, inside each checkm directory the user can find `checkm.csv` file with completeness and contamination of the re-assembled refined bins.

```
cat checkm.csv  | head -5
Bin Id,Marker lineage,# genomes,# markers,# marker sets,0,1,2,3,4,5+,Completeness,Contamination,Strain heterogeneity
bin.1,c__Deltaproteobacteria (UID3216),83,247,155,17,227,3,0,0,0,96.09,1.94,66.67
bin.2,k__Bacteria (UID2570),433,273,183,112,160,1,0,0,0,55.46,0.27,0.00
bin.3,o__Alteromonadales (UID4837),30,868,324,49,818,1,0,0,0,95.62,0.04,0.00
bin.4,f__Rhodocyclaceae (UID3972),30,540,241,132,390,18,0,0,0,77.00,2.41,11.11
```
In case the user chooses to do the filtering after **CheckM**, inside the checkm directories there is also `filtered_checkm.tsv` file where the filtered results of checkm are saved.

Each of the annotation folders contains the following directories: 

```
bin_funct_annotations  bin_translated_genes  bin_untranslated_genes  prokka_out
```
The functional annotation of the bins can be found in GFF form. As an example running: `cat Annotation_folder_bacteria/bin_funct_annotations/bin.1.gff | head -5`

```
NODE_2_length_360965	Prodigal:2.6	CDS	660	1022	.	-	0ID=EDFJOLLJ_00001;inference=ab initio prediction:Prodigal:2.6;locus_tag=EDFJOLLJ_00001;product=hypothetical protein
NODE_2_length_360965	Prodigal:2.6	CDS	1019	1306	.	-	0ID=EDFJOLLJ_00002;inference=ab initio prediction:Prodigal:2.6;locus_tag=EDFJOLLJ_00002;product=hypothetical protein
NODE_2_length_360965	Prodigal:2.6	CDS	1428	2660	.	-	0ID=EDFJOLLJ_00003;eC_number=2.6.1.83;Name=dapL;gene=dapL;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:A3DK17;locus_tag=EDFJOLLJ_00003;product=LL-diaminopimelate aminotransferase
NODE_2_length_360965	Prodigal:2.6	CDS	2816	3616	.	-	0ID=EDFJOLLJ_00004;eC_number=1.17.1.8;Name=dapB;gene=dapB;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P38103;locus_tag=EDFJOLLJ_00004;product=4-hydroxy-tetrahydrodipicolinate reductase
NODE_2_length_360965	Prodigal:2.6	CDS	3638	4510	.	-	0ID=EDFJOLLJ_00005;eC_number=4.3.3.7;Name=dapA;gene=dapA;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:O67216;locus_tag=EDFJOLLJ_00005;product=4-hydroxy-tetrahydrodipicolinate synthase

```

For more detailed explanation of the MetaWrap tool the user can study the ![metaWrap/Usage_tutorial](https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md) file. 

# Module 3: Recovery of Viral Metagenome-Assembled Genomes (wish to be finished)
Note: Make sure that all the viral tools are installed. The links for the installation can be found in the installation module of the ![README](https://github.com/mdsufz/MuDoGeR/blob/master/README.md) file.

## 3.1
By using the Module 3 the user will recover Viral Metagenome-Assembled Genomes. Also, by using Module 3 the user can estimate the quality of the viral contigs do taxonomic classification of the viruses. the taxonomy recovers on the assembly fasta file, leads to the identification and recovery of the viral genomes contained in that. The indepentent results of each tool, combined and dereplicated. Also, by using Module 3, the user can estimates the quality of the dereplicated viral contigs, do the taxonomic classification of the viruses. Furthermore, the user can choose to utilize Module 3 for determining the host of each virus. Before running the script, it is important for the user to decide about the parameters **minimum coverege(-c)** and **minimum identity (-i)** used in the dereplication. By default the minimum coverage is 70% and the minimum identity 95%. However the user is free to change the dereplication parameters depending on the aims of the metagenomic analysis or the assembled dataset. 

Run MuDoGeR viral  module:
``` 
mudoger viral module  -o /path/to/output/folder -f ~/path/to/assembly/file --g /path/to/prokaryotic_hosts/folder -c 70 -i 95 --wish
```
* The `/path/to/output/folder` indicates the path to the output directory where the output folders of the viral module will be written.
* The `/path/to/assembly/file` indicates the path to the file of the assemblies. 
* The `path/to/prokaryotic_hosts/folder` indicates the path to the directory that contains the genomes of the possible prokaryotic hosts (optional).
* The `-c` indicates that the default value of the minimum coverege in the de-replication step is 70.
* The `-i` indicates that the default value of the minimum identity in the de-replication step is 95.


In the output directory five folders are present. The `initial_recovey_folder`  contains the results from the independent recovery of each tool. The `dereplication_folder` includes the dereplication results while in the `taxonomy_folder` and `quality_folder` the user can find the results of the taxonomic classification utlizing **vContact2** tool and quality control using **CheckV** tool, respectively. In case the user chooses to use the **WiSH** tool, a fifth folder will be present, called `wish_folder`. This folder will contain the results of the host identification analysis.

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

The quality summary of the contigs can be found inside the file `quality_folder/quality_summary.csv`
Running:

```
head -5 quality_folder/quality_summary.csv 
contig_id,contig_length,genome_copies,gene_count,viral_genes,host_genes,checkv_quality,miuvig_quality,completeness,completeness_method,contamination,prophage,termini
NODE_6_length_9839_cov_5.151165,9839,1.0,16,5,0,Low-quality,Genome-fragment,16.25,AAI-based,0.0,No,No
NODE_7_length_8843_cov_4.647815,8843,1.0,10,5,0,Low-quality,Genome-fragment,15.0,AAI-based,0.0,No,No
NODE_9_length_8285_cov_3.135358,8285,1.0,14,1,0,Not-determined,Genome-fragment,NA,NA,0.0,No,No
NODE_21_length_5441_cov_6.763832,5441,1.0,8,8,0,Medium-quality,Genome-fragment,89.8,AAI-based,0.0,No,55-bp-DTR
```
In the `taxonomy_folder` the user can find  an overview of the viral clusters in `taxonomy/viral_cluster_overview.csv` file.

# Module 4: Recovery of Eukaryotic Metagenome-Assembled Genomes
Note: Make sure that all the eukaryotic tools are installed. The links for the installation can be found in the installation module of the ![README](https://github.com/mdsufz/MuDoGeR/blob/master/README.md) file.

The Eukaryotic module on the Assembly fasta file, leads to the recovery of Eukaryotic Metagenome-Assembled Genomes.

## 4.1
Running the first part of Eukaryotic module:
```
mudoger eukaryotic module 1 -f ~/path/to/assembly/file --prokarya /path/to/prokaryotic/folder -o /path/to/output/folder -1 ~/path/to/final_pure_reads_1.fastq -2 ~/path/to/final_pure_reads_2.fastq 
```

* The `/path/to/assembly/file` indicates the path to the file of the assemblies. 
* The `/path/to/output/folder` indicates the path to the output directory where the output folders of the **4.1** will be saved.
* The `/path/to/prokaryotic/folder` indicates the path to a directory where the prokaryotic assemblis will fall after the separation of eukaryotic and prokaryotic assemblies with **EukRep**. 
* The `/path/to/final_pure_reads_1.fastq` indicates the path to the file of the forward clean reads. 
* The `/path/to/final_pure_reads_2.fastq` indicates the path to the file of the reverse clean reads. 

In the output of the first step the user can find `euk_concoct_bins` folder which contains the bins produced by **CONCOCT** with size >= 2.5 Mb.

## 4.2
The second step of Eukaryotic module starts with **GeneMark-EV** and as input the user can use any of the bins produced in the previous step.

Running the second step of Eukaryotic module:

```
mudoger eukaryotic module 2 -f ~/path/to/concoct/bin/fasta/file -o /path/to/output/folder 
```
* The `/path/to/concoct/bin/fasta/file` indicates the path to the bin file.
* The `/path/to/output/folder` indicates the path to the output directory where the output folders of **4.2** will be written.

After the end of the second step the output folder contains the results from **MAKER2**, **BUSCO** and **EukCC** tools:

The results of the **EukCC** tool are located in the  `eukcc/workfiles/gmes/output/gmhmm.mod` file. The `eukcc` directory is found inside the initial output directory.

The results of the **MAKER2** tool are located in the `maker/euk-ebin.maker.output/OUTPUT.all.maker.genemark.transcripts.fasta` file which contains the names and the sequences of the annotated proteins of the predicted genes . The maker directory is found  inside the `genemark` directory. 

The results of the **BUSCO** tool are located in the file `maker/busco/full_table_fbusco.tsv`. The busco directory is found inside the `maker` directory. 





