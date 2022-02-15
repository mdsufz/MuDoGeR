
# inputs
# 1 is the taxonomy file from gtdb
# 2 is the bins(mags) folder
# 3 is the output folder where files will be dumped

taxonomy_file="$1"
bbtools_file="$2"
bins_folder="$3"
output_folder="$4"

if [ ! -d "$output_folder" ]; then mkdir "$output_folder"; else :; fi
if [ ! -d "$output_folder"/n50_good_bins ]; then mkdir "$output_folder"/n50_good_bins; else :; fi

#1 filter by n50 > 10
#while read l; do  if [ "$(echo "$l" | cut -f3 )" -gt 10 ]; then cp "$bins_folder"/"$(echo "$l" | cut -f1)"  "$output_folder"/n50_good_bins ; else :; fi; done < "$bbtools_file" 

#2 organize taxonomical clusters
#python /data/msb/tools/anisplitter/organize-bins-tax.py "$taxonomy_file" "$output_folder"/n50_good_bins "$output_folder/tax_clusters" 


#3 run fastani on taxonomical clusters
#module load Anaconda3 ; conda activate /data/msb/tools/fastani/fastani_conda
#for clu in "$output_folder/tax_clusters"/cluster*; do echo $clu; cd $clu;  for d in *; do for e in *; do  fastANI -q "$d" -r "$e" --fragLen 1500 --minFraction 0 -o fastani_output-"${d/.fa/}"-"${e/.fa/}"  ; done ; done ; cat fastani* > fastani_all.tsv    ; cd ..; done
echo "$PWD"


echo "setting up job array..."
#cd "$output_folder"/tax_clusters
#for clu in clu*; do cd $clu; for d in *fa; do for e in *fa; do echo -e ""$(realpath "$d")"\t"$(realpath "$e")"\t"$(echo "$PWD"/fastani_"${d/.fa/}"-"${e/.fa/}")""; done; done;  cd .. ; done  >   "$output_folder"/fastani_array

NUM_JOBS="$(wc -l "$output_folder"/fastani_array | cut -f1 -d' ')"
qsub -N fastani_JOBARRAY -t 1-"$NUM_JOBS" /data/msb/tools/anisplitter/strain_differentiation/jobarray_fastani/submission_script.sh "$output_folder"/fastani_array


#4 run anisplitter
#module load GCC/8.3.0  OpenMPI/3.1.4  load R
#cd "$output_folder"/tax_clusters
#for clu in cluster*; do echo $clu; mkdir $clu/ani_split_95;   /gpfs1/data/msb/tools/anisplitter/aniSplitter/R/aniSplitter.R -d "$clu"/ani_split_95 -f "$clu"/fastani_all.tsv -t "$clu"/taxonomy -a 95 ;done



