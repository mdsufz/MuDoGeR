
# Step 2:Metagenomic recovery of Prokaryotic genomes 

The help message for prokaryotic genomes recovery
```
help_message () {
        echo ""
        echo "MetaWRAP v=$VERSION"
        echo "Usage: metaWRAP [module]"
        echo ""
        echo "  Modules:"
        echo "  binning         Binning module (metabat, maxbin, or concoct)"
        echo "  bin_refinement  Refinement of archaea bins from binning module"
        echo "  bin_refinement  Refinement of bacteria bins from binning module"
        echo "  Quality_bins_a  Calculation of quality for archaea bins   
        echo "  Quality_bins_b  Calculation of quality for bacteria bins 
        echo "  Taxon_bins_a    Taxonomy for archaea bins   
        echo "  Quality_bins_b  Taxonomy for bacteria bins
        echo ""
        echo "  --help | -h     show this help message"
        echo "  --version | -v  show metaWRAP version"
        echo "  --show-config   show where the metawrap configuration files are stored"
        echo "";}

```
*  Running of the prokaryotic recovery module with the MetaWrap tool. As first input the Assembly.fasta file was used. In the second place  In the final place is necessary to give the output file for the  recovery.   

# Step 3: Pipelines for viral genomes recovery 
* In this step, the recovery of  the viral genomes from the Assembly.fasta file takes place. For the initial recovery of the viral genomes three tools are used: VirFinder, VirSorter and VIBRANT. The three tools are used simultaneously, although it is possible to be run independentely. After the initial run, the recovered data from each tool are filtered,combined and the replicates are been removed. 

* Run MuDoGeR viral recovery module:
* ```~/viral_working_directory/viral_recovery_script  -o /path/to/output/folder -f ~/path/to/assembly/file ``` 
* The generated output directory contains three fodlers, one for the initial recovered genomes of each tool and a fourth file that contains the results from the dereplication problems.



* Running of the viral recovery module with all the three recovery methods and the dereplication process (coverage 70 and id 95). The output directory is placed in the beginnng of the command before the assembly.fasta input.  
```
mudoger viral_module -o </path/to/outputdir> -f </path/to/assembly.fa>
```

* After the process is fully completed, the ouput folder contains the outputs from the initial recovery from each individual tool and the folder with the  dereplicated assembled files:

```
virfinder.tsv           virsorter_folder        vibrant_folder          extraction_and_dereplication_fasta 
``` 







