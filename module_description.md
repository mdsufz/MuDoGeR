# DETAILED DESCRIPTION OF EACH MODULE

## Viral module
The viral module utilizes three viral recovery tools(**VirFinder**, **VirFinder**, **VIBRANT**) for the identification and recovery of viral genomes from an assemble metagenomic dataset. In the beginning, each of the tools is used to recover independently the viral genomes from a given assembly dataset and locates them in seperate folders.  In the next step the recovered sequences are filtered using bash. The selection of proper VirFinder sequences is based on low q-value (q-value=< 0.01) and high lentgh (lentgh>= 1000 bp), the VirFinder chosen sequences are those classified into categories 1 and 2 and from VIBRANT the selected sequences are contained in the Assembly.phages_combined.fna file. It is important to note that in contrast with the other two tools, the VirFinder output recovery file contains only the headers of the assemblies. Because of that,with bash command, the headers of the uniq filtered sequences from each tool are extracted to a common fasta.file. In the next step, using the headers including in the common fasta.file,  the actual sequences from the assembly dataset are extracted transfered to a new fasta.file with python. Finally,the duplicated contigs are removed by a de-replication function. The dependencies of the dereplication step are  maximum coverage (-c 70) and maximum identity (-i 95). The values of the dereplication parameters can be changed depending on the users aim and datasets. The final output results, are depositet to a new output folder, contained inside the initial output directory.