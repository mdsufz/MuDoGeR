
# Step 3: Pipelines for viral genomes recovery 
Note: Make sure that all the viral tools are installed 
Running the viral module on the assembly fasta file, leads to the identification and recovery of the viral genomes contained in that. The indepentent results of each tool, combined together and the replicates are been removed. The three viral recovery tools (**VirFinder**, **VIBRANT**, **VirSorter**) are used simultaneously.  Before the run of the script is important for the user to decide about the parameters minimum coverege(-c) and minimum identity (-i) used in the dereplication. By default the minimum coverage is 70 and the minimum identity 95. However the user is free to change the dereplication parameters depending on the aims of the metagenomic analysis or the assembly dataset. 
Run MuDoGeR viral  module:
``` 
mudoger viral module  -o /path/to/output/folder -f ~/path/to/assembly/file -c 70 -i 95 
 ```

In the generated output directory four folders are appeared.  Three of them contain the results from the independent recovery of each tool while in the fourth the user can found the final de-replicated contigs. The independent recovered contigs of each tool can be usefull for other types of analysis. 

```
virfinder_folder   virsorter_folder    vibrant_folder    dereplication_folder
``` 
Inside the dereplication_folder are the final outputs of the viral module.

``` 
VIRAL_PARTICLES_95-70.clstr   VIRAL_PARTICLES-nucmer.out.coords   VIRAL_PARTICLES_95-70.fna   VIRAL_PARTICLES-cover.csv   
VIRAL_PARTICLES-nucmer.out.delta
```
     

# Step 1: Pre-Processing
## Step 1.1:

## Step 1.2:





