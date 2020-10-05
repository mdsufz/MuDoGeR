

(for the usage page)
# Raw read QC module 
Note: In the read QC module we wanted to run only the read trimming so we used the --skip-bmtagger flag.

All the samples of the input files are processing at the same time

For running the raw read QC module:

``` 
mudoger read_qc module --skip-bmtagger -1 /path/to/raw_reads_1.fasta -2 /path/to/raw_reads_2.fasta  -t ${NSLOTS:-1} -o /path/to/pure_reads/output/directory
```


The output directory of the read quality control module contains:
```
final_pure_reads_1.fastq    pre-QC_report
final_pure_reads_2.fastq    post-QC_report 
```


# Assembly Module

The reads are assembled with utilization of metaSPAdes option flag(in case of very large datasets it is bettere to use MegaHIT):
```
metawrap assembly -1 /path/to/pure_reads_1.fasta -2 path/to/pure_reads_2.fastq -m 200 -t 96 --use-metaspades -o /path/to/assembled_reads/output/directory
```
After the end of the assembly process, inside the output directory the user can find the folder `Assembly_folder`. Inside this folder is the assembly file called `final_assembly.fasta` and the QUAST assembly report html called `assembly_report.html`.


# Prokaryotic module
Note: Make sure that all the databases and programms required for the MetaWrap run are downloaded

The run of the prokaryotic module leads to the recovery of prokaryotic genomes from the assembly dataset by using the MetaWrap tool. The script of the prokaryotic module combines the algorithms of every MetaWrap module. The prokaryotic module starts with the binning of the co-assembled reads using the binning module 

Run the prokayotic module with using MetaWrap:
``` 
mudoger prokaryotic module -o /path/to/metawrap/output/directory -f ~/path/to/assembly/file ---metawrap -1 ~/path/to/final_pure_reads_1.fastq -2 -1 ~/path/to/final_pure_reads_2.fastq 

```
In the final output folder the user can find:
