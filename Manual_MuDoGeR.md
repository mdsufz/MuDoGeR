
# Step 3: Pipelines for viral genomes recovery 
Note: Make sure that all the viral tools are installed. 

Running the viral module on the assembly fasta file, leads to the identification and recovery of the viral genomes contained in that. The indepentent results of each tool, combined together and the replicates are been removed. The three viral recovery tools (**VirFinder**, **VIBRANT**, **VirSorter**) are used simultaneously but there is also the choice for the user to run each of them independently or even to skip one of them.  Before the run of the script is important for the user to decide about the parameters minimum coverege(-c) and minimum identity (-i) used in the dereplication. By default the minimum coverage is 70 and the minimum identity 95. However the user is free to change the dereplication parameters depending on the aims of the metagenomic analysis or the assembly dataset. 

Run MuDoGeR viral  module with all three tools and the de-replication function:
``` 
mudoger viral module  -o /path/to/output/folder -f ~/path/to/assembly/file --vifinder --virsorter --vibrant --dereplication -c 70 -i 95 
 ```

In the final output directory four folders are appeared. Three of them contain the results from the independent recovery of each tool while in the fourth the user can found the final de-replicated contigs and the . The independent recovered contigs of each tool can be usefull for other types of analysis. 

```
virfinder_folder   virsorter_folder    vibrant_folder    dereplication_folder
``` 
Inside the `dereplication_folder` are only the final output files of the viral module.

``` 
VIRAL_PARTICLES_95-70.clstr   VIRAL_PARTICLES-nucmer.out.coords   VIRAL_PARTICLES_95-70.fna   VIRAL_PARTICLES-cover.csv   
VIRAL_PARTICLES-nucmer.out.delta
``` 
The  `VIRAL_PARTICLES_95-70.clstr` file contains the header and the length of the clusters. For example `head -5 VIRAL_PARTICLES_95-70.clstr`:

```
>Cluster_0	NODE_53_length_4168_cov_5.199611	4168 
>Cluster_1	NODE_58_length_3975_cov_7.914796	3975
>Cluster_2	NODE_67_length_3721_cov_2.320513	3721
>Cluster_3	NODE_105_length_3251_cov_4.997497	3251
>Cluster_4	NODE_108_length_3205_cov_6.333333	3205

```
The sequences of the clusters can be found in the `VIRAL_PARTICLES_95-70.fna` file.









