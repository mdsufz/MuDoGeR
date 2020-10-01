

# Step 3: Pipelines for viral genomes recovery 
* In this step, the recovery of  the viral genomes from the Assembly.fasta file takes place. For the initial recovery of the viral genomes three tools are used: VirFinder, VirSorter and VIBRANT. The three tools are used simultaneously, although it is possible to be run independentely. After the initial run, the recovered data from each tool are filtered,combined and the replicates are been removed. 
* Run MuDoGeR viral recovery module:
* ```~/viral_working_directory/viral_recovery_script  -o /path/to/output/folder -f ~/path/to/assembly/file  -c 70 -i 95 ``` 
* Before the run of the script is important for the user to decide about the parameters minimum coverege(-c) and minimum identity (-i) used in the dereplication. 
By default minimum coverage is 70 and the minimum identity 95. However the user is free to change the dereplication parameters depending on the aims of the analysis

* The generated output directory contains three fodlers, one for the initial recovered genomes of each tool and a fourth file that contains the results from the dereplication problems.



* Running of the viral recovery module with all the three recovery methods and the dereplication process (coverage 70 and id 95). The output directory is placed in the beginnng of the command before the assembly.fasta input.  
```
mudoger viral_module -o </path/to/outputdir> -f </path/to/assembly.fa>
```

* After the process is fully completed, the ouput folder contains the outputs from the initial recovery from each individual tool and the folder with the  dereplicated assembled files:

```
virfinder.tsv           virsorter_folder        vibrant_folder          extraction_and_dereplication_fasta 
``` 







