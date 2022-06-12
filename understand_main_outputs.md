# Main Outputs

![ScreenShot](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/mudoger_outputs.jpg)

The Multi-Domain Genome Recovery v1.0 (MuDoGeR v1.0) framework is a tool developed to help users to recover Prokaryotic Metagenome-Assembled Genomes (pMAGs), Uncultivated Viral Genomes (UViGs), and Eukaryotic Metagenome-Assembled Bins (eMABs) from whole-genome sequence (WGS) samples simultaneously. The MuDoGeR v1.0 framework act as a wrapper of several tools. It was designed to be an easy-to-use tool that orchestrates the application of the used bioinformatics tools and outputs ready-to-use comprehensive files. You should be able to run 1 simple command for each of the five modules inside MuDoGeR.
Therefore, when using MuDoGeR you will be able to find the relevant outputs for all used tools. Please refer to the Manual to understand where to find the tools-specific outputs and information on how they were created. In the Manual, you can also find the reference to the tools themselves, where you can read their specific documentation.
 However, to make the Genome recovery analysis easier, MuDoGeR also parses and summaries the resulted files into ready-to-use comprehensive files. Following, you find the description and location of the final files created by MuDoGeR.

## Module 1 – Pre-processing and Assembly

From module 1, the 3 main outputs you will find are:

* The quality-controlled reads located at ```output_folder/sample_name/qc/final_pure_reads_1/2.fastq```. The files are in regular .fastq format;
* The RAM prediction for [**metaSPAdes**](https://genome.cshlp.org/content/27/5/824) assembly, located at ```output_folder/sample_name/khmer/final_prediction.tsv```. This is a tab-separated-values (tsv) file with only two columns: ```dataset   maxvmem_M```. They indicate the file used to predict the memory (dataset), and the memory requirements in megabytes (MB) (maxvmem_M)
* The assembled sequence in regular .fasta file located at ```output_folder/sample_name/assembly/final_assembly.fasta```

The folder structure from the results of Module 1 should be as follows:

```
sample_name
     ├── assembly
     │   ├── final_assembly.fasta
     │   └── …
     ├── khmer
     │   ├── final_prediction.tsv
     │   └── …
     └── qc
         ├── final_pure_reads_1.fastq
         ├── final_pure_reads_2.fastq
         └── …
 ```
 
## Module 2 – pMAGs recovery final files

From the recovery of pMAGs, MuDoGeR parses the outputs from the tools mentioned in the Manual and outputs the following comprehensive results to the ```sample_name/prokaryotes/final_outputs/``` folder:

* The recovered prokaryotic bins are located in the folder ```all_bins_seq/```. Those files are named as “sample_name”-bin.”number”.fa. The numbering simply indicates different sequences. They are all regular fasta files.
* The ```allbins_metrics_summary.tsv``` file. This file contains a summary of all information recovered from all prokaryotic bins. It is a regular .tsv file and its columns are: ``` OTU     completeness    contamination   str.heterogeneity       taxonomy        genome_size     #scaffolds      largest_scaff   N50     N90     prokka_know	prokka_unknown```. They are, respectively: the bin name, completeness by [**CheckM**](https://genome.cshlp.org/content/25/7/1043), contamination by [**CheckM**](https://genome.cshlp.org/content/25/7/1043), strain heterogeneity by [**CheckM**](https://genome.cshlp.org/content/25/7/1043), taxonomy by [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182), bin size, number of scaffolds, size of the largest scaffold, N50 value, N90 value, number of known annotated [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) genes, number of unknown identified [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) genes.
*  The ```bins_genes_prokka_summary/``` folder containing the gene annotation summary performed by [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) for each of the recovered bins. The tsv file contains the following columns: ``` locus_tag       ftype   length_bp       gene    EC_number       COG     product```. You can also check the [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) documentation by clicking on the tool's name.
*  The ```bins_metrics_summary/ qual_bins_checkm_summary.tsv``` file. This file contains the summary of the quality estimation procedure performed by [**CheckM**](https://genome.cshlp.org/content/25/7/1043). You can check the official [**CheckM**](https://genome.cshlp.org/content/25/7/1043) documentation by clicking on the tool's name.
*  ``` bins_metrics_summary/ taxa_bins_gtdbtk_summary.tsv ``` file. This file contains the summary of the final taxonomical classification procedure performed by [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182). You can check the official [**GTDB-tk**](https://academic.oup.com/bioinformatics/article/36/6/1925/5626182) documentation by clicking on the tool's name.
* The ```mags_results_summary.tsv``` file. This file is subset from the ```allbins_metrics_summary.tsv``` file containing only the bins classified as MAGs (completeness – 5*contamination > 50). 
* The ```only_mags_seq/``` folder containing a copy of the fasta file from the bins classified as MAGs.

```
sample_name/prokaryotes/final_outputs
              ├── all_bins_seq
              │    ├── sample_name-bin.0.fa
              │    ├── sample_name-bin.1.fa
              │    ├── sample_name-bin.2.fa
              │    └── …
              ├── allbins_metrics_summary.tsv
              ├── bins_genes_prokka_summary
              │    ├── sample_name-bin.0_genes_prokka.tsv
              │    ├── sample_name-bin.1_genes_prokka.tsv
              │    ├── sample_name-bin.2_genes_prokka.tsv
              │    └── …
              ├── bins_metrics_summary
              │    ├── qual_bins_checkm_summary.tsv
              │    └── taxa_bins_gtdbtk_summary.tsv
              ├── mags_results_summary.tsv
              └── only_mags_seq
                   ├── sample_name-bin.2.fa
                   ├── sample_name-bin.4.fa
                   └── …
````

## Module 3 – UViGs recovery final files

From the recovery of UViGs, MuDoGeR parses the outputs from the tools mentioned in the Manual and outputs the following comprehensive results to the ```sample_name/viruses/final_outputs/``` folder:

* The ``` only_uvigs_seq/``` folder contains the viral sequences classified as UViGs following the standard defined by [Roux, S., et al.(2019)](https://www.nature.com/articles/nbt.4306) and [Nayfach, S., et al. (2021)](https://www.nature.com/articles/s41587-020-00774-7). The files are regular fasta files and are named as “sample_name”_putative_viral_contig-“number”.fa
* The ``` putative_vir_contigs_summary.tsv``` file. This file is a summary of all retrieved information from all the putative_vir_contigs recovered during the complete MuDoGeR pipeline. It is a regular tsv file, and its columns are: ``` uvig    original_contig uvig_length     provirus        proviral_length gene_count      viral_genes     host_genes     checkv_quality  miuvig_quality  completeness    completeness_method     contamination   kmer_freq     warnings putative_host   likelihood```. The columns mean, respectively, the given name to the uvig, the original contig name from the uvig found in the .fasta file, the uvig length in bp, provirus classification, proviral length (if any), number of genes found, number of viral genes identified, number of identified host genes, quality classification by [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7), miuvig quality classification, completeness, completeness method used, khmer frequency, any warnings identified, name of the prokaryotic host identified by [**WiSH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377), and probability of that prokaryotic being the host.
* The ```putative_vir_seq_metrics_summary/host_vir_pair_wish_summary.tsv``` file. This file contains the summary of the Viral-Host pair estimation performed by [**WiSH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377). You can check the official [**WiSH**](https://academic.oup.com/bioinformatics/article/33/19/3113/3964377) documentation by clicking on the tool's name.
* The ```putative_vir_seq_metrics_summary/quality_checkv_summary.tsv``` file. This file contains the summary of the quality estimation performed by [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7). You can check the official [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7) documentation by clicking on the tool's name.
* The ```putative_vir_seq_metrics_summary/taxa_estimation_vir_vcontact2_summary.csv``` file. This file contains a parsed summary of the taxonomical classification performed by vcontact. To check how MuDoGeR parses the [**Vcontact2**](https://www.nature.com/articles/s41587-019-0100-8) output, please, check MuDoGeR’s Manual. You can check the official vcontact documentation by clicking on the tool's name. The columns of this file are: ``` uvig, original_contig, Domain, Order, Family ```. The columns indicate, respectively, the uvig name according to the file names created by MuDoGeR, the original contig name from the uvig found with the .fasta file, the domain taxa level, the order taxa level, and the family taxa level.
* The ```putative_vir_seq_metrics_summary/vir_contigs_annotation_vibrant_summary.tsv``` file. This file contains the summary of the viral gene annotation performed by [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0) on all putative viral contigs. You can check the official [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0) documentation by clicking on the tool's name.
* The ```putative_vir_seq_metrics_summary/vir_contigs_genbank_annotation_vibrant.tsv``` file. This file contains the Genbank summary annotation performed by [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0). You can check the official [**VIBRANT**](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00867-0) documentation by clicking on the tool's name.
* The ``` putative_viral_contigs/``` folder containing regular fasta sequences files from all recovered putative viral sequences recovered in MuDoGeR.
* The ``` uvigs_high_quality_summary.tsv``` file. This file is a regular tsv and it is a subset from the ``` putative_vir_contigs_summary.tsv``` containing only complete and High-quality sequences as classified by [**CheckV**](https://www.nature.com/articles/s41587-020-00774-7).
* The ```viral_contigs_seq_names.csv``` file. This is a regular CSV file containing a direct mapping from the MuDoGeR putative viral sequences names to the original contig names found within the fasta files.

The folder structure from the results of Module 3 should be as follows:

```
sample_name/viruses/final_outputs
           ├── only_uvigs_seq
           │    ├── sample_name_putative_viral_contig-15.fa
           │    ├── sample_name_putative_viral_contig-16.fa
           │    └── ...
           ├── putative_vir_contigs_summary.tsv
           ├── putative_vir_seq_metrics_summary
           │    ├── host_vir_pair_wish_summary.tsv
           │    ├── quality_checkv_summary.tsv
           │    ├── taxa_estimation_vir_vcontact2_summary.csv
           │    ├── vir_contigs_annotation_vibrant_summary.tsv
           │    └── vir_contigs_genbank_annotation_vibrant.tsv
           ├── putative_viral_contigs
           │    ├── sample_name_putative_viral_contig-0.fa
           │    ├── sample_name_putative_viral_contig-1.fa
           │    └── …
           ├── uvigs_high_quality_summary.tsv
           └── viral_contigs_seq_names.csv
```


## Module 4 – eMABs recovery final files




## Module 5 – pMAGs/UViGs/eMABs coverage and relative abundance tables

Initially, on step 5.a, MuDoGeR selects the representative pMAGs using the gOTUpick software. If you are interested in knowing the MAGs belonging to the same OTU group, you can find this information in the ```mapping_results/gOTUpick_results/``` folder. There you will find the results files:

* The ```bestbins.txt``` file. This a two-column file containing the name of the representative bin in the first column and the group it represents in the second. Bins from unique taxonomies are classified as the best bins from that taxonomical group.
* The ```final_groups_outputs.csv``` file. This is a regular CSV file with 3 columns. The first indicates the name of all grouped bins, the second indicates the group that bin was assigned, and the third indicates with a “\*” if that bin the representative from its group. 

The results from step 5.a are structured as follows:
```
mapping_results/gOTUpick_results
            ├── final_output
            ├── bestbins.txt
            └── final_groups_output.csv
```

Step 5.b performs the mapping and calculation of the coverage and relative abundance from the recovered representative eMABs, pMAGs, and UViGs. After running step 5.b MuDoGeR creates the ```mapping_results/euk_mabs_mappings``` ```mapping_results/pmags_otu_mapping/```, and ```mapping_results/uvigs_mapping/``` for storing the results from the eMABs, pMAGs, and UViGs respectively. Inside those folders you should find the ```map_final_tables_reduced/``` folder if you use the MuDoGeR with the ```--reduced``` flag, and the ```map_final_tables_complete/``` if you use the ```--complete``` flag. Inside those folders you will find the final output files created by MuDoGeR as follows: 

* The ```map_reduced(complete)_absolute_n_hits_list.tsv``` file. This is a regular three-column tsv file containing in the first column the eMABs, pMAGs, and UViGs, the second is the sample used for the mapping, and the third is the absolute number of mapped reads. This file is created if you use the ```--reduced (--complete)``` flag.
* The ``` map_reduced(complete)_absolute_n_hits_table.tsv``` file. This file is created based on the ```map_reduced(complete)_absolute_n_hits_list.tsv``` file. This is a regular tsv file with samples used for mapping as columns and target eMABs, pMAGs, and UViGs as rows. 
* The ```map_reduced(complete)_coverage_list.tsv``` file. This is a regular three-column tsv file containing in the first column the target eMABs, pMAGs, and UViGs, the second is the sample used for the mapping, and the third is the coverage calculated for the eMABs, pMAGs, and UViGs. This file is created if you use the ```--reduced (--complete)``` flag.
* The ``` map_reduced(complete)_coverage_table.tsv``` file. This file is created based on the ```map_reduced(complete)_coverage_list.tsv``` file. This is a regular tsv file with samples used for mapping as columns and target eMABs, pMAGs, and UViGs as rows. 
* The ```map_reduced(complete)_relative_abundance_list.tsv``` file. This is a regular three-column tsv file containing in the first column the target eMABs, pMAGs, and UViGs, the second is the sample used for the mapping, and the third is the relative abundance of the eMABs, pMAGs, and UViGs on the used samples. This file is created if you use the ```--reduced (--complete)``` flag.
* The ```map_reduced(complete)_relative_abundance_table.tsv``` file. This file is created based on the ```map_reduced(complete)_relative_abundance_list.tsv``` file. This is a regular tsv file with samples used for mapping as columns and target eMABs, pMAGs, and UViGs as rows. 

The results from step 5.b are structured as follows:
```
mapping_results/(euk_mabs_mapping)(pmags_otu_mapping)(uvigs_mapping)/
                    ├── map_final_tables_complete
                    │	    ├── map_complete_absolute_n_hits_list.tsv
                    │	    ├── map_complete_absolute_n_hits_table.tsv
                    │	    ├── map_complete_coverage_list.tsv
                    │	    ├── map_complete_coverage_table.tsv
                    │	    ├── map_complete_relative_abundance_list.tsv
                    │	    └── map_complete_relative_abundance_table.tsv
                    └── map_final_tables_reduced
                        ├── map_reduced_absolute_n_hits_list.tsv
                        ├── map_reduced_absolute_n_hits_table.tsv
                        ├── map_reduced_coverage_list.tsv
                        ├── map_reduced_coverage_table.tsv
                        ├── map_reduced_relative_abundance_list.tsv
                        └── map_reduced_relative_abundance_table.tsv
```

Step 5.c performs the genes' relative abundance calculation from the samples' assembly. After using MuDoGeR module 5 using the ```--gene``` flag you should have the folders ```mapping_results/assembly_gene_map/```. Inside that folder you should find the following most relevant output files: 

* The ```functional_annotation/sample_name/``` folder. There you will find the [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517) gene annotation performed on the ```final_assembly.fasta``` file from the ```sample_name``` library.
* The ```genelength/sample_name.genelength``` file. This is a two-column file with the gene id in the first and the gene length in the second from the genes annotated by [**Prokka**](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517).
* The ``` map_absolute_count/sample_name.count``` file. This is a two-column file with the gene id in the first and the absolute number of reads mapped to that gene.
* The ``` map_coverage_norm/sample_name.cov``` file. This is a two-column file with the gene id in the first and the gene coverage on the sample_name where that gene was annotated. The coverage is calculated as: (average read length)\*(number of hits)/(gene length)
* The ``` map_tpm_norm/sample_name.tpm``` file. This is a two-column file with the gene id in the first and the gene TPM calculation on the sample_name where that gene was annotated. The TPM calculation is done as: (reads mapped)/(average read length)\*10^6/(gene length)\*Sum\[(reads mapped)\*(average read length)/(gene length)\].

The folder structure from the results of 5.c should be as follows:

```
mapping_results/ assembly_gene_map/
                     ├── functional_annotation
                     │   	└── sample_name
                     ├── genelength
                     │   	└── sample_name.genelength
                     ├── map_absolute_count
                     │   	└── sample_name.count
                     ├── map_coverage_norm
                     │   	└── sample_name.cov
                     └── map_tpm_norm
                          └── sample_name.tpm
```



# Using the tools independently

A complete metagenome recovery pipeline can require several bioinformatics tools with complex installation and usage procedures. The MuDoGeR pipeline makes it easy to install and run a complete metagenome-assembled genome recovery.  Consequently, the MuDoGeR pipeline installs several tools that have multiple dependencies. More often than not, those tools require conflicting dependencies. To solve this problem, MuDoGeR will use miniconda to create multiple environments that are automatically orchestrated during the pipeline's functioning.
Therefore, after installing MuDoGeR the user will have a working conda environment for each of the integrated tools. If the user wants to customize the use of any tool, you can use MuDoGeR to configure your machine and activate its environments independently. 
The MuDoGeR installation procedure uses the variable ```$CONDA_PREFIX``` to detect your base environment, so make sure you have conda properly installed on your system. During installation, MuDoGeR saves several relevant system paths from your system on the file called ```config.sh```. In addition, MuDoGeR creates a mudoger_env and configures several new environments within the mudoger_env directory. The environments' names are called ```toolname_env```. 
The MuDoGeR configures and installs are the following environments named according to their main tool:

```
bbtools_env		brat_env		busco_env		checkv_env
eukcc_env		eukrep_env		extract_env		genemarker_env
gtdbtk_env		htseq_env		khmer_env		maker2_env
metawrap_env		otupick_env		prokka_env
stampede_clustergenomes_env		vcontact2_env		vibrant_env
virfinder_env		virsorter2_env		wish_env
```

Essentially, the user can activate a tool-specific environment by running:
```console
$ conda activate path/to/your/miniconda/envs/mudoger_env/dependencies/conda/envs/toolname_env
```

However, since the user’s conda path is also saved by MuDoGeR during installation, the user can also activate a tool-specific environment by doing the following:

```console
$ conda activate mudoger_env
$ config_path="$(which config.sh)"
$ database="${config_path/config/database}"
$ source $config_path
$ source $database
$ conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/toolname_env
````

