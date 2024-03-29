#!/bin/bash
# PROKARYOTES MODULES FILESYSTEM STRUCTURE
#   BINNING
#     INITIAL-BINNING
#     REFINEMENT-BACTERIA
#     REFINEMENT-ARCHAEA
#     BIN-FINAL-SET (REDUNDANCY REMOVAL)
#   METRICS
#     CHECKM
#     GTDBtk
#     PROKKA
#     STATS (GENOMIC)
#       N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more...


help_message () {
	echo ""
	echo "Usage: bash -i MuDoGeR/src/modules/mudoger-module-2.sh -1 reads_1.fastq -2 reads_2.fastq -o output_dir -a assembly.fasta"
	echo "Options:"
	echo ""
	echo "	-1 STR			forward fastq reads path"
	echo "	-2 STR			reverse fastq reads path"
	echo "	-a STR			assembly path"
	echo "	-o STR			output directory path"
	echo "	-m INT			given RAM in GB (default=50G)"
	echo "	-t INT			number of threads/cores (default=1)"
	echo "	-h --help		print this message"
	echo "";}

#DEFINE DEFAULT PARAMETERS

# the option memory was desabled as we are using the memory predicted by script 1.2. We need to develop another way if it does not work with the amount of memory predicted. 

libname_folder=$(pwd) 		#output path for the downloaded sequences
memory=50			#given Memory to the Assembly process in GB
cores=1 			#number of threads that is going to be used

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


echo 'running'

conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path


# 1 INITIAL BINNING USING METAWRAP (CONCOCT, METABAT2, MAXBIN2)
if [ -d  "$libname_folder"/prokaryotes/binning/initial-binning/concoct_bins ]; 		# if one of the outputs is already there, do not run
then echo "-> Binning already done. Please check here: "$libname_folder"/prokaryotes/binning/initial-binning"
else
echo "-> Running initial binning"
mkdir -p "$libname_folder"/prokaryotes/binning
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-1_initial-binning.sh "$assembly"            \
                                      "$forward_library"                                  \
                                      "$reverse_library"                                   \
                                      "$libname_folder"/prokaryotes/binning/initial-binning \
                                      "$cores"                                               \  #
				      "$memory"
fi

# 2 BIN REFINEMENT USING METAWRAP FOR BACTERIA (50,10)
if [ -f  "$libname_folder"/prokaryotes/binning/refinement-bac/metawrap_50_10_bins.stats ];
then echo "-> Bacterial refinement is done. Please check here: "$libname_folder"/prokaryotes/binning/refinement-bac"
else
echo "-> Running bacterial refinement"
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-2_bin-ref-bacteria.sh "$libname_folder"/prokaryotes/binning/refinement-bac  \
                                       "$cores"                                                         \
                                       "$assembly"                                                       \
                                       "$libname_folder"/prokaryotes/binning/initial-binning/concoct_bins \
                                       "$libname_folder"/prokaryotes/binning/initial-binning/maxbin2_bins  \
                                       "$libname_folder"/prokaryotes/binning/initial-binning/metabat2_bins  \
                                       "$memory"       
fi
              

# 3 BIN REFINEMENT USING METAWRAP FOR ARCHAEA (40,30)
if [ -f  "$libname_folder"/prokaryotes/binning/refinement-arc/metawrap_40_30_bins.stats ];
then echo "-> Archaeal refinement is done. Please check here: "$libname_folder"/prokaryotes/binning/refinement-arc"
else
echo "-> Running archaeal refinement"
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-3_bin-ref-archea.sh "$libname_folder"/prokaryotes/binning/refinement-arc \
                                     "$cores"                                                        \
                                     "$assembly"                                                       \
                                     "$libname_folder"/prokaryotes/binning/initial-binning/concoct_bins \
                                     "$libname_folder"/prokaryotes/binning/initial-binning/maxbin2_bins   \
                                     "$libname_folder"/prokaryotes/binning/initial-binning/metabat2_bins   \
				     "$memory" 
fi


# 4 BIN REDUNANCY REMOVAL
if [ -d  "$libname_folder"/prokaryotes/binning/unique_bins ];
then echo "-> Bin dereplication is done. Please check: "$libname_folder"/prokaryotes/binning/unique_bins"
else
echo '-> Run bin redundancy removal.'
mkdir -p "$libname_folder"/prokaryotes/binning/unique_bins
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-4_bin-dereplication.sh  "$libname_folder" 
fi


# 5 GTDBtk taxonomy assignment
mkdir -p "$libname_folder"/prokaryotes/metrics
if [ -f "$libname_folder"/prokaryotes/metrics/GTDBtk_taxonomy/gtdbtk_result.tsv ];
then echo "-> Bin taxonomy assignment is done. Please check: "$libname_folder"/prokaryotes/metrics/GTDBtk_taxonomy"
else
echo "-> Run bin taxonomy"
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-5_bin-taxonomy.sh "$libname_folder"/prokaryotes "$cores"
fi


# 6 CheckM quality control
if [ -f "$libname_folder"/prokaryotes/metrics/checkm_qc/outputcheckm.tsv ];
then echo "-> MAGs quality control is done. Please check: "$libname_folder"/prokaryotes/metrics/checkm_qc/outputcheckm.tsv"
else 
echo "-> Run mags checkm"
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-6_bin-QC.sh "$libname_folder"/prokaryotes/ "$cores"
fi

# 7 PROKKA Annotation
if [ -d "$libname_folder"/prokaryotes/metrics/prokka ] ; 
then echo "-> Rapid prokka annotation is done. Please check: "$libname_folder"/prokaryotes/metrics/prokka"
else 
echo "Run Prokka annotation"
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-7_prokka.sh  "$libname_folder"/prokaryotes/ "$cores"
fi


# 8 Metrics: N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more...

if [ -f "$libname_folder"/prokaryotes/metrics/genome_statistics/genome_metrics.tsv ] && [ "$libname_folder"/prokaryotes/metrics/genome_statistics/prok_genomes_stats.tsv ];
then echo "-> MAGs statistics is done. Please check: "$libname_folder"/prokaryotes/metrics/genome_statistics/prok_genomes_stats.tsv"
else 
echo "-> Run MAGs STATS"
bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-8_genomics-metrics.sh  "$libname_folder"/prokaryotes/ "$cores"
fi


