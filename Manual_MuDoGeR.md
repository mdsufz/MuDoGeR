

# Step 3: Pipelines for viral genomes recovery 
* Note: Make sure that all the viral tools are installed 
* Running the viral module on the assembly fasta file, leads to the identification and recovery of the viral genomes contained in that. The indepentent results of each tool, combined together and the replicates are been removed. The three viral recovery tools (**VirFinder**, **VIBRANT**, **VirSorter**) are used simultaneously, but it is possible for each tool to be run seperately and the user is possible to chose not to run one of the tools, depending on the aim of the analysis and the size of the assembly dataset.

* Run MuDoGeR viral  module:


* ```~/viral_working_directory/viral_recovery_script  -o /path/to/output/folder -f ~/path/to/assembly/file   ``` 
* Before the run of the script is important for the user to decide about the parameters minimum coverege(-c) and minimum identity (-i) used in the dereplication. 
By default minimum coverage is 70 and the minimum identity 95. However the user is free to change the dereplication parameters depending on the aims of the analysis

* The generated output directory contains three fodlers, one for the initial recovered genomes of each tool and a fourth file that contains the results from the dereplication problems.



* Running of the viral recovery module with all the three recovery methods and the dereplication process (coverage 70 and id 95). The output directory is placed in the beginnng of the command before the assembly.fasta input.  
* ```mudoger viral_module -o </path/to/outputdir> -f </path/to/assembly.fa> ```

* After the process is fully completed, the ouput folder contains the outputs from the initial recovery from each individual tool and the folder with the  dereplicated assembled files:

```
virfinder.tsv           virsorter_folder        vibrant_folder          extraction_and_dereplication_fasta 
``` 







