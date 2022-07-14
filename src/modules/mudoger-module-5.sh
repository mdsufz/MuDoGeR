#!/bin/bash

#BRAT: a tool to calculate relative abundance of fasta files in fastq files
#MERGE RESULTS and run BRAT
#Running bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/$module_script -o $output_folder -t $num_cores $brat_type;
#The user should have initial data to use BRAT: output from GTDB-Tk, output from CheckM, output from BBMap, the recovered bins/FASTA in fasta format, and Quality Controlled Pair End reads.

help_message() {
		echo ""
		echo "Example usage: BRAT.sh  -f path/to/fasta_files_folder -r path/to/reads_fastq_folder -o path/to/outputdir --reduced --coverage"
		echo ""
		echo "You must select one BRAT type: --reduced  or --complete BRAT"
		echo ""
		echo "You must select at least one output option: --absolute-values and/or --coverage and/or --relative-abundance"
		echo ""
		echo "Available options:"
		echo ""
		echo "-h, --help"
		echo "	display this help page"
		echo ""
		echo "--meta <value>"
		echo "	absolute path to metadata file as explained in MuDoGeR github documentation"
		echo "-f <value>"
		echo "	absolute path to a directory where fasta files are are stored (.fa)"
		echo ""
		echo "-r <value>"
		echo "	absolute path to a directory where reads in fastq format are stored (.fastq)"
		echo ""
		echo "-o <value>"
		echo "	absolute path to a directory where the user wants to store the output"
		echo ""
		echo ""
		echo "-t <value>"
		echo "	number of threads/cores [default: 1]"
		echo ""
		echo "BRAT types:"
		echo ""
		echo "--reduced"	
		echo "	run BRAT of the fasta files only to the origin’s read. For prokaryoten, it will be the representative OTU to the reads of all other OTU inside the same species’ group"
		echo ""
		echo "--complete"	
		echo "	run BRAT of the fasta files to all reads."
		echo ""
		echo "--genes"	
		echo "	Annotate genes and run BRAT on the assembled contigs."
		echo ""
		echo "Output types:"
		echo ""
		echo "--absolute-values"	
		echo "	calculate absolute number of reads mapped."
		echo ""
		echo "--coverage"	
		echo "	calculate coverage."
		echo ""
		echo "--relative-abundance"	
		echo "	calculate relative abundance."
		echo "";
				}

#Default parameters
#"WORKDIR" is the path to the output directory.
#"FASTA" is the path to the directory where FASTA/bins are found. 
#"READS" is the path to the directory where samples are found. 

WORKDIR=None;
THREAD=1;
reduced=true; #default
complete=false;
absolute=true; #default
coverage=false;
relative=false;
genes=false;

OPTS='getopt -o f:r:o:t --long reduced:,complete:,absolute-values:,coverage:,relative-abundance:,help -- "$@" '

if [ $? -ne 0 ]; then echo "Warning: Something unexpected happened" help_message; exit 1; fi

while true 
do
	case $1 in
		-o) WORKDIR="$2"; shift 2;;
		-t) THREAD="$2"; shift 2;;
		--meta) metadata_table=$2; shift 2;;
		--reduced) reduced=true; shift 1;; 
		--complete) complete=true;reduced=false; shift 1;;
		--genes) genes=true; reduced=false; complete=false; shift 1;; 
		--absolute-values) absolute=true; shift 1;; 
		--coverage) coverage=true; shift 1;; 
		--relative-abundance) relative=true; shift 1;; 
		-h | --help) help_message; exit 1; shift 1;;
		*) break;;
	esac
done

#Checking if output path are entered correctly
if [ "$WORKDIR" = "None" ] || [ "$metadata_table" = "None" ] ; then 
	echo "Non-optional parameter for BRAT type was not entered. Please select output path -o path/to/output/folder and metadata table"
	help_message; exit 1
fi

#Checking if all non-optional parameters are entered correctly
if [ "$reduced" = "false" ] && [ "$complete" = "false" ] && [ "$genes" = "false" ]; then 
	echo "Non-optional parameter for BRAT type was not entered. Please select --reduced or --complete."
	help_message; exit 1
fi

#Checking if all non-optional parameters are entered correctly
if [ "$reduced" = "true" ] && [ "$complete" = "true" ]; then 
	echo "Non-optional parameter for BRAT type was not entered. Please only select one option, --reduced or --complete."
	help_message; exit 1
fi

#Checking if all non-optional parameters are entered correctly
if [ "$absolute" = "false" ] && [ "$coverage" = "false" ] && [ "$relative" = "false" ]; then 
	echo "Non-optional parameters for output type were not entered. Please only select one option, --absolute-values ,--coverage or--relative-abundance."
	help_message; exit 1
fi

#Checking if output directory already exists
if [ ! -d $WORKDIR/mapping_results ]; then mkdir -p $WORKDIR/mapping_results; fi


echo -e "\nTHE PIPELINE STARTED\n"


#Module 5 - BRATs for MuDoGer

conda activate mudoger_env
config_path="$(which config.sh)"
source $config_path

#	OTUpick                         submodule 5-1
#	Bin Mapping process             submodule 5-2
#	Gene mapping on assembly	submodule 5-3

echo -e "\n MAP CALCULATION STARTED"

#   OTU picking      submodule 5-1
if [ -f $WORKDIR/mapping_results/gOTUpick_results/final_output/bestbins.txt ]; then
	echo -e "\nGenome OTU pick already done. Check '$WORKDIR'/mapping_results/gOTUpick_results/final_output/bestbins.txt";
else
	echo -e "\n OTU picking STARTED"
	bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-5-1_prokOTUpicking.sh $WORKDIR $metadata_table $THREAD
	echo -e "\n OTU picking DONE"
fi
	

if [ "$reduced" = "true" ] || [ "$complete" = "true" ]; then

	#     Perform BRAT      submodule 5-2
	
	echo -e "\n BRAT STARTED"
	bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-5-2_BRAT.sh $WORKDIR \
										$metadata_table \
										$THREAD \
										$absolute \
										$coverage \
										$relative \
										$reduced \
										$complete
	echo -e "\n BRAT DONE"


	
	
elif [ "$genes" = "true" ]; then

	#Gene mapping on assembly	submodule 5-3
	echo -e "\n Gene mapping STARTED"
	bash -i $MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-5-3_genemap.sh $WORKDIR \
										$metadata_table \
										$THREAD \
										$absolute \
										$coverage \
										$relative
										
	echo -e "\n Gene mapping DONE"



fi


echo -e "\nThat's the end of BRAT. Enjoy your results. Goodbye!"
