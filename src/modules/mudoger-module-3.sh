#!/bin/bash
# VIRUSES MODULES FILESYSTEM STRUCTURE
#   VIRAL INVESTIGATION
#     VIRSORTER               submodule 3-1
#     VIRFINDER               submodule 3-1
#     VIBRANT                 submodule 3-1
#     DEREPLICATION           submodule 3-1
#   METRICS
#     VCONTACT2 (TAXONOMY)     submodule 3-2
#     WISH (HOST PREDICTION)   submodule 3-3
#     VCHECK (QUALITY)         submodule 3-4
#     STATS (GENOMIC)
#       N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more... submodule 3-5


help_message () {
	echo ""
	echo "Usage: bash -i MuDoGeR/src/modules/mudoger-module-2.sh -1 reads_1.fastq -2 reads_2.fastq -o output_dir -a assembly.fasta"
	echo "Options:"
	echo ""
	echo "	-1 STR			forward fastq reads path"
	echo "	-2 STR			reverse fastq reads path"
	echo "	-a STR			assembly path"
	echo "	-o STR			output directory path"
	echo "	-m INT			given Memory to the Assembly process in GB (default=10)"
	echo "	-t INT			number of threads/cores (default=1)"
	echo "	-h --help		print this message"
	echo "";}

# loop through input params
while true; do
	case "$1" in
		-1) forward_library=$2; shift 2;;
		-2) reverse_library=$2; shift 2;;
		-a) assembly=$2; shift 2;;
		-o) libname_folder=$2; shift 2;;
		-t) cores=$2; shift 2;;
		-m) memory=$2; shift 2;;
		-h | --help) help_message; exit 1; shift 1;;
		--) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done

conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path

#assembly=$1
#libname_folder=$2              # master output folder to be created and defined by user
#cores=10

# 1 viral investigation
#mkdir -p "$libname_folder"/viruses/investigation
#bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-3-1_viral-investigation.sh "$assembly" "$libname_folder"/viruses/investigation $cores

if [ -f  "$libname_folder"/viruses/investigation/vibrant_filtered_data.txt ] && [ -f  "$libname_folder"/viruses/investigation/virfinder_filtered_data.txt ] && [ -f  "$libname_folder"/viruses/investigation/virsorter2_filtered_data.txt ]; 		# if the outputs is already there, do not run
then echo "-> Viral investigation already done. Please check here: "$libname_folder"/viruses/investigation/"
else
echo "-> Running Viral investigation"
mkdir -p "$libname_folder"/viruses/investigation
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-3-1_viral-investigation.sh "$assembly" \
"$libname_folder"/viruses/investigation \
$cores

fi


# 2 viral vcontact2 (taxonomy)
#bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-3-2_vcontact2.sh "$libname_folder"/viruses/investigation/dereplication/uvigs_95-70.fna "$libname_folder"/viruses/taxonomy "$cores"

if [ -f  "$libname_folder"/viruses/taxonomy/vcontact-output/genome_by_genome_overview.csv ]; 		# if one of the outputs is already there, do not run
then echo "-> Viral taxonomy already done. Please check here: "$libname_folder"/viruses/taxonomy/vcontact-output/genome_by_genome_overview.csv"
else
echo "-> Running Viral taxonomy"

bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-3-2_vcontact2.sh "$libname_folder"/viruses/investigation/dereplication/uvigs_95-70.fna \
"$libname_folder"/viruses/taxonomy \
"$cores"

fi

# 3 host prediction

#bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-3-3_host-prediction.sh "$libname_folder"/viruses/investigation/dereplication/uvigs_95-70.fna  "$libname_folder"/prokaryotes/binning/unique_bins "$libname_folder"/viruses/host_prediction

if [ -f  "$libname_folder"/viruses/host_prediction/output_results/prediction.list ]; 		# if one of the outputs is already there, do not run
then echo "-> Running Viral Host Prediction already done. Please check here: "$libname_folder"/viruses/host_prediction/output_results/prediction.list"
else
echo "-> Running Viral Host Prediction"

bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-3-3_host-prediction.sh "$libname_folder"/viruses/investigation/dereplication/uvigs_95-70.fna \
"$libname_folder"/prokaryotes/binning/unique_bins \
"$libname_folder"/viruses/host_prediction

fi

# 4 vcheck

#bash -i MuDoGeR/src/scripts/mudoger-module-3-4_vcheck.sh "$libname_folder"/viruses/investigation/dereplication/uvigs_95-70.fna "$libname_folder"/viruses/vcheck_quality "$cores"
if [ -d  "$libname_folder"/viruses/vcheck_quality/quality_summary.tsv ]; 		# if one of the outputs is already there, do not run
then echo "-> Viral quality check already done. Please check here: "$libname_folder"/viruses/vcheck_quality/quality_summary.tsv"
else
echo "-> Running Viral quality check"

bash -i MuDoGeR/src/scripts/mudoger-module-3-4_vcheck.sh "$libname_folder"/viruses/investigation/dereplication/uvigs_95-70.fna \
"$libname_folder"/viruses/vcheck_quality \
"$cores"

fi

# 5 uvigs metrics
#bash -i MuDoGeR/src/scripts/mudoger-module-3-5_uvigs-metrics.sh


