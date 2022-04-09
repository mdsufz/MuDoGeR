#!/bin/bash

#gOTUpick: a tool to identify strains from the genomes of closely related species

#The user should have an initial data to use gOTUpick: output from GTDB-Tk, output from CheckM, output from BBMap and the recovered bins/MAGs in fasta format.

help_message() {
		echo ""
		echo "Example usage: gOTUpick.sh --bb-input path/to/BBMap-output --checkm-input path/to/CheckM-output --gtdb-input path/to/gtdb-output 
				      -m path/to/mags -o path/to/outputdir" 
		echo ""
		echo ""
		echo "Available options:"
		echo ""
		echo "-h, --help"
		echo "	display this help page"
		echo ""
		echo "--bb-input <value>"
		echo "	absolute path to output file from BBMap"
		echo ""
		echo "--checkm-input <value>"
		echo "	absolute path to output file from CheckM"
		echo ""
		echo "--gtdb-input <value>"
		echo "	absolute path to output file from GTDB-Tk"
		echo ""
		echo "-m <value>"
		echo "	absolute path to a directory where mags/bins in fasta format are stored"
		echo ""
		echo "-o <value>"
		echo "	absolute path to a directory where the user wants to store the output"
		echo ""
		echo "--fastANI-fragLen <value>"
		echo "	fragment lenght for fastANI [default: 1500]"
		echo ""
		echo "--fastANI-minFraction <value>"	
		echo "	minFraction for fastANI [default: 0]"
		echo ""
		echo "--fastANI-thread <value>"
		echo "	thread count for paralellization in fastANI [default: 2]"
		echo ""
		echo "-s <value>"
		echo "	set seed for reproducibility [default: 2020]"
		echo ""
		echo "--a1 <value>" 
		echo "	ANI threshold for splitting groups into clusters [default: 95]"
		echo ""
		echo "--a2 <value>"
		echo "	ANI threshold for splitting clusters into subclusters (i.e. second splitting) [default: 99]"
		echo ""
		echo "-b <value>"	
		echo "	number of iterations for bootstraping [default: 1000]"
		echo ""
		echo "--no-prefilter"	
		echo "	do not filter out low-quality mags/bins"
		echo "";
				}

#Default parameters
#"WORKDIR" is path to output directpry. "BB" is path to BBMap input. "CHECKM" is path to CheckM input. 
#"GTDB" is path to GTDB-Tk input. "MAGS" is path to the directory where mags/bins are found. "FRAGLEN" is the fragment length in fastANI.
#"MINFRAC" is the minimum fraction of genome in fastANI. "THREAD" is the thread number for parellization in fastANI.
#"SEED" is an arbitrary number set for reproducibility of the run. "ANI_1" is the ANI threshold used in first clustering.
#"ANI_2" is the ANI threshold used in second clustering. "BOOTSTRAP" is the number of iterations for bootstraping.
#"prefilter" defines whether low quality mags/bins will be filtered out. 

#Manage envs and path
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database

conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env

WORKDIR=None; BB=None; CHECKM=None; GTDB=None; MAGS=None; FRAGLEN=1500; MINFRAC=0; THREAD=2; SEED=2020; ANI_1=95; ANI_2=99; BOOTSTRAP=1000; prefilter=true

OPTS='getopt -o m:o:s:b:h --long bb-input:,checkm-input:,gtdb-input:,fastANI-fragLen:,fastANI-minFraction:,fastANI-thread:,a1:,a2:,
no-prefilter,help -- "$@"'

if [ $? -ne 0 ]; then echo "Warning: Something unexpected happened" help_message; exit 1; fi

while true 
do
	case $1 in
		--bb-input) BB="$2"; shift 2;;
		--checkm-input) CHECKM="$2"; shift 2;;
		--gtdb-input) GTDB="$2"; shift 2;;
		-m) MAGS="$2"; shift 2;;
		-o) WORKDIR="$2"; shift 2;;
		--fastANI-fragLen) FRAGLEN="$2"; shift 2;;
		--fastANI-minFraction) MINFRAC="$2"; shift 2;;
		--fastANI-thread) THREAD="$2"; shift 2;;
		-s) SEED="$2"; shift 2;;
		--a1) ANI_1="$2"; shift 2;;
		--a2) ANI_2="$2"; shift 2;;
		-b) BOOTSTRAP="$2"; shift 2;;
		--no-prefilter) prefilter=false; shift 1;; 
		-h | --help) help_message; exit 1; shift 1;;
		*) break;;
	esac
done

#Checking if all non-optional parameters are entered correctly
if [ "$WORKDIR" = "None" ] || [ "$BB" = "None" ] ||  [ "$CHECKM" = "None" ] || [ "$GTDB" = "None" ] || [ "$MAGS" = "None" ]; then 
	echo "Non-optional parameters --bb-input and/or --checkm-input and/or --gtdb-input were not entered"
	help_message; exit 1
fi

#Checking if output directory already exists
if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi

#Checking if BBMap file with a content inside exists
if [ ! -s $BB ]; then 
	echo "Warning: Input file from BBMap does not exist or has no content inside. Exiting..."; exit 1; 
fi

#Checking if CheckM file with a content inside exists
if [ ! -s $CHECKM ]; then 
	echo "Warning: Input file from CheckM does not exist or has no content inside. Exiting..."; exit 1; 
fi

#Checking if GTDB file exists and has content inside
if [ ! -s $GTDB ]; then 
	echo "Warning: file from GTDB does not exist or has no content inside. Exiting..."; exit 1; 
fi

#Checking if directory for mags exists and has files inside
if [ ! -d $MAGS ]; then
	echo "Warning: The directory for mags/bins does not exist . Exiting..."; exit 1;
else 
	cd $MAGS
	if [ "$(ls *fa | wc -l)" -eq 0 ]; then
		echo "Warning: There is no mag/bin inside mags/bins directory. Exiting..."; exit 1;
	fi
fi

cd $WORKDIR

echo -e "\nTHE PIPELINE STARTS\n"

########DATA PREPARATION AND GROUPING BY TAXONOMY########

#Sorting input files according to the mag/bin name
head -n+1 $BB > .header_bbtools.tsv
tail -n+2 $BB | sort -k20,20 > .sorted_content_bbtools.tsv
cat .header_bbtools.tsv .sorted_content_bbtools.tsv > sorted_bbtools.tsv
mv sorted_bbtools.tsv $BB
rm -f .header_bbtools.tsv .sorted_content_bbtools.tsv

head -n+1 $CHECKM > .header_output-checkm.tsv
tail -n+2 $CHECKM | sort -k1,1 > .sorted_content_output-checkm.tsv
cat .header_output-checkm.tsv .sorted_content_output-checkm.tsv > sorted_output-checkm.tsv
mv sorted_output-checkm.tsv $CHECKM
rm -f .header_output-checkm.tsv .sorted_content_output-checkm.tsv

head -n+1 $GTDB > .header_gtdb_taxonomy.tsv
tail -n+2 $GTDB | sort -k1,1 > .sorted_content_gtdb_taxonomy.tsv
cat .header_gtdb_taxonomy.tsv .sorted_content_gtdb_taxonomy.tsv > sorted_gtdb_taxonomy.tsv
mv sorted_gtdb_taxonomy.tsv $GTDB
rm -f .header_gtdb_taxonomy.tsv .sorted_content_gtdb_taxonomy.tsv



#Creating the directory for taxonomy groups. Later, the bins/MAGs with the same GTDG-Tk taxonomy will be stored in a single sub-folder in the tax_groups.
mkdir -p $WORKDIR/tax_groups


#Pre-filtering of low-quality recovered genomes. The condition is that N50 is greater than 10000 and quality score is equal or greater than 50.
#The following block of code filters out low-quality genomes, only if the user did not disable it by choosing the option "no-prefilter".

if [ "$prefilter" = "true" ]; then
	paste <(tail -n+2 $BB) <(tail -n+2 $CHECKM) | awk 'BEGIN {FS="\t"}; $7 >= 10000 && ($32 - (5*$33)) >= 50 {print $20}' | rev | cut -d '/' -f1 | sed 's/af.//1' | rev > goodqual.ids
	
	#Creating the a new GTDB file which has the taxonomy information only for good-quality bins/MAGs.
	while read l; 
	do 
		grep -P "$l\t" $GTDB; 
	done < goodqual.ids > gtdb_taxonomy_qc.tsv
	
	TAXFILE="gtdb_taxonomy_qc.tsv"; 

else
	
	tail -n+2 $GTDB > allgtdb.tsv
	
	TAXFILE="allgtdb.tsv"
fi

#echo -e "\nPlease enter your anaconda path (ex: /home/user_name/anaconda3):" #needed to find where auxiliary scripts are located 
#read conda_path
#echo -e "\nYour anaconda path: $conda_path" 

echo -e "\nStarting first step: Grouping by taxonomy"

#The Python script organize-bins-tax.py groups the good-quality bins/MAGs based on the taxonomy assigned by GTDB-Tk.
#python3 $conda_path/envs/gOTUpick/bin/organize-bins-tax.py $TAXFILE $MAGS $WORKDIR/tax_groups
python3 "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env/bin/organize-bins-tax.py $TAXFILE $MAGS $WORKDIR/tax_groups
cd -
cd $WORKDIR/tax_groups

if [ "$(ls | wc -l)" -eq 0 ]; then
	echo "Warning: Grouping by taxonomy did not work succesfully. Exiting..."; exit 1
else 
	for i in $WORKDIR/tax_groups/*
	do 
		cd $i
		if [ "$(ls *fa | wc -l)" -eq 0 ]; then
			echo "Warning: Grouping by taxonomy did not work succesfully. Exiting..."; exit 1
		fi
	done
fi
	
echo -e "\nGrouping by taxonomy is finished. Starting next step: Calculation of ANI distances\n"
cd -
########CALCULATION OF ANI DISTANCES########

cd $WORKDIR

#Creating a folder within the working directory to store the fastANI output.
mkdir -p $WORKDIR/ANI_distances

#Creating files containing paths to the bins/MAGs. The pathnames will be the input for fastANI.
for i in $WORKDIR/tax_groups/gr*
do 
	cid=$(echo $i | rev | cut -f1 -d '/' | rev); mkdir -p $WORKDIR/ANI_distances/$cid; realpath $i/*.fa > $WORKDIR/ANI_distances/$cid/paths_$cid.txt
done

echo "fastANI is starting..."

#Running fastANI to calculate ANI distances 
for i in ANI_distances/gr*
do 
	fastANI --ql $i/paths* --rl $i/paths* -o $i/fastani-out-1500-0.txt --fragLen $FRAGLEN --minFraction $MINFRAC -t $THREAD
	if [ ! -s $i/fastani-out-1500-0.txt ]; then echo "Warning: Something went wrong with fastANI. Exiting..."; exit 1; fi
done

#Copying GTDB-taxonomy files for each group into the folders in ANI_distances directory
cd -
cd ANI_distances

for i in gr*
do 
	cp ../tax_groups/$i/group_taxonomy_gtdb $i/. 
done

echo -e "\nCalculation of ANI distances is finished. Starting next step: Clustering\n"
cd -
########CLUSTERING MAGS/BINS INSIDE EACH GROUP USING CALCULATED ANI DISTANCES########

#Running aniSplitter.R with threshold ANI=95 (this is the default value). 
#aniSplitter clusters bins/mags in each taxonomy group based on ANI distances using hierarchical agglomerative clustering.

cd $WORKDIR

for i in $WORKDIR/ANI_distances/gr*
do 
	#$conda_path/envs/gOTUpick/bin/aniSplitter.R -t group_taxonomy_gtdb -d $i -f fastani-out-1500-0.txt -s $SEED -a $ANI_1 -i $BOOTSTRAP 2> $i/log.txt
	Rscript "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env/bin/aniSplitter.R -t group_taxonomy_gtdb -d $i -f fastani-out-1500-0.txt -s $SEED -a $ANI_1 -i $BOOTSTRAP 2> $i/log.txt
	if [ ! -s $i/cluster_summary.tsv ]; then echo "Warning: Something went wrong with aniSplitter. Exiting..."; exit 1; fi
done

#aniSplitter was designed to do generate one output per folder. Remember that you might have more than one folder (i.e. tax group).
#To summarize the results of multiple outputs we employ the python script summarize-anisplitter-results.py.
cd -
mkdir -p $WORKDIR/results
cd $WORKDIR/ANI_distances

for i in gr*
do
	#python3 $conda_path/envs/gOTUpick/bin/summarize-anisplitter-results.py $i $ANI_1
	python3 "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env/bin/summarize-anisplitter-results.py $i $ANI_1
done > $WORKDIR/results/groups_ANI$ANI_1.tsv

if [ ! -s $WORKDIR/results/groups_ANI$ANI_1.tsv ]; then echo "Warning: Something went wrong with summarizing aniSplitter results. Exiting..."; exit 1; fi

sed -i '1i original_group\tbin\tbootstrap\tnew_group\tani\tbootstrap_status' $WORKDIR/results/groups_ANI$ANI_1.tsv
cd -
cd $WORKDIR/

#Combining all ANI results (those obtained from fastANI) in a single file
cat $WORKDIR/ANI_distances/gr*/fastani-out* > $WORKDIR/all_fastani-out-1500-0.txt

#Creating a new directory called ANI_OTU95
mkdir -p $WORKDIR/ANI_OTU$ANI_1
cd -
cd $WORKDIR/ANI_OTU$ANI_1

#Creating folders in the directory ANI_OTU95 for the clusters obtained after the running anisplitter with the threshold ANI=95 
cut -f4 $WORKDIR/results/groups_ANI$ANI_1.tsv | tail -n+2 | sort | uniq | while read l; do mkdir $l; done;

#Copying GTDB-taxonomies to the folders for respective clusters 
for i in gr*
do 
	orig=$(echo $i | cut -f1-2 -d '-'); echo $orig; cp $WORKDIR/ANI_distances/$orig/group_taxonomy_gtdb $i
done

#Copying ANI distances to the folders for clusters: each ANI distance file will contain pairwise ANI distances only for the bins/mags in the respective cluster.
for i in gr*
do 
	for j in $(grep "$i\s" <(tail -n+2 $WORKDIR/results/groups_ANI$ANI_1.tsv) | cut -f2 )
	do 
		grep $j $WORKDIR/all_fastani-out-1500-0.txt | grep -vf <(grep -v "$i\s" <(tail -n+2 $WORKDIR/results/groups_ANI$ANI_1.tsv) | cut -f2) >> $i/fastani-out-1500-0.txt; echo $i
	done
done

#Removing duplicates in the ANI distance files.
for i in gr*
do 
	sort $i/fastani-out-1500-0.txt | uniq > $i/tmp; mv -f $i/tmp $i/fastani-out-1500-0.txt
done

#Identifying new clusters with only one sequence and moving these clusters into a folder called unique
#Note that the folder unique is created only if there is a cluster with one sequence

for i in gr*/fastani-out-1500-0.txt
do 
	wc -l $i
done | grep "^1 " | rev | cut -f2 -d '/' | rev | cut -f2 -d " " | grep -v "OUT"| while read l
do 
	mkdir -p unique
	mv $l $WORKDIR/ANI_OTU$ANI_1/unique/
done

if  [ ! -d "unique" ]; then
	echo "No cluster with only one sequencce"
fi

#The above command line moves the clusters with only one sequence into the folder "unique". One exception is when the clusters with only 
#one sequence has bad bootstrap values. In this case, these clusters are not moved into to the folder "unique", but, into the folder "out" 
#by using the following command lines, together with the other clusters having bad bootstrap values (if there are any).

#If there are clusters that have bad bootstrap values (bootstrap < 75), create a separate folder for them and move them to this folder
if [ "$(ls -d *OUT | wc -l)" -gt 0 ]; then
	mkdir -p $WORKDIR/ANI_OTU$ANI_1/out
	mv $WORKDIR/ANI_OTU$ANI_1/gr*OUT $WORKDIR/ANI_OTU$ANI_1/out/;

else 
	echo "No cluster with bad bootstrap values"
fi

#Running ANI splitter with the threshold ANI=99 (default) on clusters in the ANI_OTU95 directory
for i in gr*
do 
	#$conda_path/envs/gOTUpick/bin/aniSplitter.R -t group_taxonomy_gtdb -d $i -f fastani-out-1500-0.txt -s $SEED -a $ANI_2 -i $BOOTSTRAP 2> $i/log.txt
	Rscript "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env/bin/aniSplitter.R -t group_taxonomy_gtdb -d $i -f fastani-out-1500-0.txt -s $SEED -a $ANI_2 -i $BOOTSTRAP 2> $i/log.txt
	if [ ! -s $i/cluster_summary.tsv ]; then echo "Warning: Something went wrong with aniSplitter. Exiting..."; exit 1; fi
done
cd -
#Running summarize-anisplitter-results.py script using a threshold of 99 and storing the results in the result folder.
cd $WORKDIR/ANI_OTU$ANI_1

for i in gr*
do 
	#python3 $conda_path/envs/gOTUpick/bin/summarize-anisplitter-results.py $i $ANI_2
	python3 "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env/bin/summarize-anisplitter-results.py $i $ANI_2
done > $WORKDIR/results/groups_ANI$ANI_2.tsv

if [ ! -s $WORKDIR/results/groups_ANI$ANI_2.tsv ]; then echo "Warning: Something went wrong with summarizing aniSplitter results. Exiting..."; exit 1; fi

sed -i '1i original_group\tbin\tbootstrap\tnew_group\tani\tbootstrap_status' $WORKDIR/results/groups_ANI$ANI_2.tsv
cd -
#Comparing groups_95 and groups_99 we need to know which clusters were further divided at 99 and which were not,then create a definitive file with the new "groups"

#Figuring out which bins did not get new clusters at 99

#First change your working directory to "results"
cd $WORKDIR/results

#Identifying the clusters that were not further subdivided using a threshold of 99
#Identifying the clusters that were subdivided using a threshold of 99
#Combining all the clusters generated into a single file

cut -f4 $WORKDIR/results/groups_ANI$ANI_2.tsv | tail -n+2 | sort | uniq | cut -f1,2,3 -d '-' | sort | uniq -c | sed -re 's/\s+//' | grep  ^1 | cut -f2 -d ' ' > $WORKDIR/results/not_further_$ANI_1.txt
grep -vf $WORKDIR/results/not_further_$ANI_1.txt groups_ANI$ANI_2.tsv | cut -f4 | tail -n+2 | cut -f1,2,3 -d '-' | sort | uniq > $WORKDIR/results/further_$ANI_1.txt
cat <(grep -vf further_$ANI_1.txt groups_ANI$ANI_1.tsv) <(grep -vf not_further_$ANI_1.txt groups_ANI$ANI_2.tsv | tail -n+2) > $WORKDIR/results/final_groups.tsv


#Generating files to extract metadata. Those files are based on the "checkm",  "bbmap" and "gtdb" outputs

#Selecting metadata for checkm output ordered

while read l
do 
	grep -P "${l/.fa/}\t" $CHECKM
done < <(tail -n+2 final_groups.tsv | cut -f2) > metadata_checkm_sorted.tsv

if [ ! -s $WORKDIR/results/metadata_checkm_sorted.tsv ]; then echo "Warning: Something went wrong while creating CheckM metadata file. Exiting..."; exit 1; fi

#Calculating and adding quality to last column
perl -ni -e 'chomp; @f=split("\t"); $qual=$f[11]-(5*$f[12]); print("$_"."\t"."$qual"."\n")' metadata_checkm_sorted.tsv

#Selecting metadata for bbmap output ordered
while read l
do 
	grep -P "${l/.fa/}$" $BB
done < <(tail -n+2 final_groups.tsv | cut -f2) > metadata_bbtools_sorted.tsv

if [ ! -s $WORKDIR/results/metadata_bbtools_sorted.tsv ]; then echo "Warning: Something went wrong while creating BBMap metadata file. Exiting..."; exit 1; fi

#Selecting metadata for gtdb-tk output ordered
while read l
do 
	grep -P "${l/.fa/}\t" $GTDB
done < <(tail -n+2 final_groups.tsv | cut -f2) > metadata_gtdb_taxonomy_sorted.tsv

if [ ! -s $WORKDIR/results/metadata_gtdb_taxonomy_sorted.tsv ]; then echo "Warning: Something went wrong while creating GTDB metadata file. Exiting..."; exit 1; fi


#Combining files to final groups
paste <(tail -n+2 final_groups.tsv | cut -f2) <(tail -n+2 final_groups.tsv | cut -f1,3,4,5,6) <(cut -f12,13,14,15 metadata_checkm_sorted.tsv) <(cut -f2,7 metadata_bbtools_sorted.tsv) <(cut -f2 metadata_gtdb_taxonomy_sorted.tsv) > final_groups_qual.tsv

if [ ! -s $WORKDIR/results/final_groups_qual.tsv ]; then echo "Warning: Something went wrong while combining metadata. Exiting..."; exit 1; fi

sed -i '1i bin\toriginal_group\tbootstrap\tnew_group\tani\tbootstrap_status\tcompleteness\tcontamination\tstrain_heterogeneity\tquality_score\tcontig_number\tn50\ttaxonomy' final_groups_qual.tsv


echo -e "\nClustering is finished. Starting next step: Selecting representative mags/bins\n"
cd -
########SELECTING REPRESENTATIVE MAGS/BINS######

#In this part, the best bin per cluster will be selected. To do that, the highest quality score (completeness - 5 x contamination) is checked. 
#If there is a tie, the lowest number of contigs, then the highest N50 and lastly, the lowest strain heterogeneity is checked.

cd $WORKDIR/results

#Running the script "pick_rep.R" to chose the best bins; the output will be a file containing best bins (bestbins.tsv)

#Rscript $conda_path/envs/gOTUpick/bin/pick_rep.R $WORKDIR/results/final_groups_qual.tsv
Rscript "$MUDOGER_DEPENDENCIES_ENVS_PATH"/otupick_env/bin/pick_rep.R $WORKDIR/results/final_groups_qual.tsv

if [ ! -s $WORKDIR/results/bestbins_metadata.tsv ]; then echo "Something went wrong while picking representatives per cluster. Exiting..."; exit 1; fi 

#In the beginning, while grouping bins/mags by taxonomy, there might have been a case where only 1 bin has a given taxonomy. 
#If this is the case, these bins are stored in the folder unique_tax. 
#Therefore, in the end we need to remember to add the unique tax that were not even included in the ANI analysis.
#They will be representative bins of their own taxonomy because they are unique.
#Note that the folder unique_tax was only created if there are bins with unique taxonomy.
cd -


mkdir -p $WORKDIR/final_output #Creating the folder to store final output

if [ -d "$WORKDIR/tax_groups/unique_tax" ]; then 

	#Changing working directory to where your unique taxonomic bins are located
	cd $WORKDIR/tax_groups/unique_tax 
	#Adding the unique taxonomic bins to the bestbins.tsv file
	cat <(cat $WORKDIR/results/bestbins_metadata.tsv | cut -f1,4 | tail -n+2) <(for i in *;do echo -e "$i\tunique taxonomy"; done) > $WORKDIR/final_output/bestbins.txt 
	#Adding header to the file
	sed -i '1i bin\tcluster_id' $WORKDIR/final_output/bestbins.txt

else
	cd $WORKDIR
	cat $WORKDIR/bestbins_metadata.tsv | cut -f1,4 | tail -n+2 > $WORKDIR/final_output/bestbins.txt
	sed -i '1i bin\tcluster_id' $WORKDIR/final_output/bestbins.txt
fi

if [ ! -s $WORKDIR/final_output/bestbins.txt ]; then echo "Something went wrong while creating final output file. Exiting..."; exit 1; fi 
cd -
echo -e "\nThat's the end of gOTUpick. Enjoy your clusters and representative bins. Goodbye!"
