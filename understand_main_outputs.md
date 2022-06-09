# Main Outputs

![ScreenShot](https://github.com/mdsufz/MuDoGeR/blob/master/flowcharts/mudoger_outputs.jpg)

The Multi-Domain Genome Recovery v1.0 (MuDoGeR v1.0) framework is a tool developed to help users to recover Prokaryotic Metagenome-Assembled Genomes (pMAGs), Uncultivated Viral Genomes (UViGs), and Eukaryotic Metagenome-Assembled Bins (eMABs) from whole-genome sequence (WGS) samples simultaneously. The MuDoGeR v1.0 framework act as a wrapper of several tools. It was designed to be an easy-to-use tool that orchestrates the application of the used bioinformatics tools and outputs ready-to-use comprehensive files. You should be able to run 1 simple command for each of the five modules inside MuDoGeR.
Therefore, when using MuDoGeR you will be able to find the relevant outputs for all used tools. Please refer to the Manual to understand where to find the tools-specific outputs and information on how they were created. In the Manual, you can also find the reference to the tools themselves, where you can read their specific documentation.
 However, to make the Genome recovery analysis easier, MuDoGeR also parses and summaries the resulted files into ready-to-use comprehensive files. Following, you find the description and location of the final files created by MuDoGeR.

## Module 1 – Pre-processing and Assembly

From module 1, the 3 main outputs you will find are:

* The quality-controlled reads located at ```output_folder/sample_name/qc/final_pure_reads_1/2.fastq```. The files are in regular .fastq format;
* The RAM prediction for metaSPAdes assembly, located at ```output_folder/sample_name/khmer/final_prediction.tsv```. This is a tab-separated-values (tsv) file with only two columns: ```dataset   maxvmem_M```. They indicate the file used to predict the memory (dataset), and the memory requirements in megabytes (MB) (maxvmem_M)
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

```
sample_name/prokaryotes/final_outputs
├── all_bins_seq
│   ├── A-bin.0.fa
│   ├── A-bin.1.fa
│   ├── A-bin.2.fa
│   └── …
├── allbins_metrics_summary.tsv
├── bins_genes_prokka_summary
│   ├── A-bin.0_genes_prokka.tsv
│   ├── A-bin.1_genes_prokka.tsv
│   ├── A-bin.2_genes_prokka.tsv
│   └── …
├── bins_metrics_summary
│   ├── qual_bins_checkm_summary.tsv
│   └── taxa_bins_gtdbtk_summary.tsv
├── mags_results_summary.tsv
└── only_mags_seq
    ├── A-bin.2.fa
    ├── A-bin.4.fa
    └── …
```

## Module 3 – UViGs recovery final files

```
sample_name/viruses/final_outputs
├── only_uvigs_seq
│   ├── A_putative_viral_contig-15.fa
│   ├── A_putative_viral_contig-16.fa
│   └── ...
├── putative_vir_contigs_summary.tsv
├── putative_vir_seq_metrics_summary
│   ├── host_vir_pair_wish_summary.tsv
│   ├── quality_checkv_summary.tsv
│   ├── taxa_estimation_vir_vcontact2_summary.csv
│   ├── vir_contigs_annotation_vibrant_summary.tsv
│   └── vir_contigs_genbank_annotation_vibrant.tsv
├── putative_viral_contigs
│   ├── A_putative_viral_contig-0.fa
│   ├── A_putative_viral_contig-1.fa
│   └── …
├── uvigs_high_quality_summary.tsv
└── viral_contigs_seq_names.csv
```

## Module 4 – eMABs recovery final files



## Module 5 – pMAGs/UViGs/eMABs coverage and relative abundance tables
```
mapping_results/gOTUpick_results
  ├── final_output
  	├── bestbins.txt
 └── final_groups_output.csv
```

```
mapping_results/ euk_mabs_mapping(pmags_otu_mapping)(uvigs_mapping)/
├── map_final_tables_complete
│	   ├── map_complete_absolute_n_hits_list.tsv
│	   ├── map_complete_absolute_n_hits_table.tsv
│	   ├── map_complete_coverage_list.tsv
│	   ├── map_complete_coverage_table.tsv
│	   ├── map_complete_relative_abundance_list.tsv
│	   └── map_complete_relative_abundance_table.tsv
└── map_final_tables_reduced
   ├── map_reduced_absolute_n_hits_list.tsv
   ├── map_reduced_absolute_n_hits_table.tsv
   ├── map_reduced_coverage_list.tsv
   ├── map_reduced_coverage_table.tsv
   ├── map_reduced_relative_abundance_list.tsv
   └── map_reduced_relative_abundance_table.tsv
```

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
