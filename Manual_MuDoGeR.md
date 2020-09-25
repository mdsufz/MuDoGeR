
# Step 3: Pipelines for viral genomes recovery

The help message of this step:
```
help_message () {
        echo ""
        echo "Usage: Viral assembly and dereplication [options] -f Assembly.fasta -o output_dir"
        echo "Options:"
        echo ""
        echo "  -f STR          Assembly.fasta"
        echo "  -o STR          output directory"          
        echo ""
        echo "  -virfinder     Recovery of viral data with VirFinder"
        echo "  -virsorter     Recovery of viral data with VirSorter"
        echo "  -vibrant       Recovery of viral data with VIBRANT"
        echo "  -dereplication Removal of replicate sequences"
        echo "";}

```
* Running of the viral module with all the three recovery methods and the dereplication process (coverage 70 and id 95). The output directory is placed in the beginnng of the command before the assembly.fasta inpu. In the final place is necessary to give the output file for the VirFinder recovery.  
* ```mudoger viral_module -o </path/to/outputdir> -f </path/to/assembly.fa> </path/to/output/folder/file/virfinder.tsv>```

* After the process is fully completed, the ouput folder contains the outputs from the initial recovery from each individual tool and the folder with the  dereplicated assembled files:

```
virfinder.tsv  virsorter_folder vibrant_folder dereplication_outputs 
`` 







