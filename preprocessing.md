

(for the usage page)
# Raw read QC module 
Note: In the read QC module we wanted to run only the read trimming so we used the --skip-bmtagger flag.

All the samples of the input files are processing at the same time

For running the raw read QC module:

``` 
mudoger read_qc module --skip-bmtagger -1 /path/to/raw_reads_1.fasta -2 /path/to/raw_reads_2.fasta  -t ${NSLOTS:-1} -o /path/to/pure_reads/output/directory
```


In the output directory 


# ASSEMBLY Module


# Prokaryotic module
Note: Make sure that all the databases and programms required for the MetaWrap run are downloaded

The run of the prokaryotic module leads to the recovery of prokaryotic genomes from the assembly dataset by using the MetaWrap tool. The script of the prokaryotic module combines the algorithms of every MetaWrap modules. The prokaryotic module start with the binning of the co-assemly process starts with binnning of the co-assembly using the binning module 

Run the prokayotic module with using MetaWrap:
``` 
mudoger prokaryotic module -o /path/to/metawrap/output/directory -f ~/path/to/assembly/file ---metawrap -1 ~/path/to/final_pure_reads_1.fastq -2 -1 ~/path/to/final_pure_reads_2.fastq 

``` 

