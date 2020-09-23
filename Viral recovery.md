
# Step 3:Pipelines for viral genomes recovery(VirFinder, VirSorter, Vibrant)


* Run the pipelines for viral genomes recovery module ``` ~/viral_working_directory/viral_recovery_script -f ~/path/to/assembly/file -o /path/to/output/folder -o /path/to/output/folder/file/virsorter.tsv``` 

* Note: To achieve revovery and derepeplication process the modules mentioned in the scirpt should be loaded. In different case the script for the viral ricavery will not work.

### 3.1 Viral genomes recovery **VirFinder**, **VirSorter**, **VIBRANT**,:

* Create the output directory:
```mkdir $1
output_viral="$1" 
``` 
* VirFinder

```
Rscript /data/msb/tools/virfinder/virfinder_script.r ${output_viral} $2 $3
```
* VirSorter 

```
output_virsorter="$output_viral/virsorter"

wrapper_phage_contigs_sorter_iPlant.pl -f $2 --wdir $output_virsorter --ncpu ${NSLOTS:-1}
#VIBRANT 
(not ready yet, needs correction)
```

* VIBRANT (still have something to be corrected, not ready yet)


### 3.2 Filtering of the results, combination and removal of repeated sequences   

* Filtering 
```#VirFinder
cat $output_viral/virfinder.tsv | awk -F'\t' '{ if ( $4 <= 0.01) print }' | awk -F'_' '{ if ( $4 >= 1000) print  }' | cut -f2 | sed "s/\"//g" > $output_viral/vir4/virf

#VirSorter
cat $virsorter_filt_inp/Predicted_viral_sequences/VIRSorter_cat-{1..2}*fasta | grep ">" | sed "s/>VIRSorter_//g"  | sed "s/-cat_2//g" | sed "s/-cat_1//g" | sed 's/\(.*\)_/\1./' > $output_viral/vir4/virs_filt

#VIBRANT (not sure yet because of the correction needed in the VIBRAN. But as I know what result gives, I know what will be the files and how to filtering it. #They will be  some changes for sure but the  main idea is this)

vibrant_filt_inp="$output_viral/vibrant_file"
#command (needs correction later)
cat $vibrant_filt_inp/*phages*combined*fna | grep ">" | sed "s/_fragment_1//g;s/>//g" > $output_viral/vibr_filt 
```

* Combination of results and removal the repeates sequences

``` cd $output_viral/vir4
cat * virs_filt vibr_filt virf | sort | uniq > COMBINED_VIRAL_PARTICLES_FOR_EXTRACTION
```

### 3.3 Extraction of the  


