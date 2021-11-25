#!/bin/bash

mudoger_lettering(){
echo -e "\n"
echo -e "\t███    ███ ██    ██ ██████   ██████   ██████  ███████ ██████  "
echo -e "\t████  ████ ██    ██ ██   ██ ██    ██ ██       ██      ██   ██ "
echo -e "\t██ ████ ██ ██    ██ ██   ██ ██    ██ ██   ███ █████   ██████  "
echo -e "\t██  ██  ██ ██    ██ ██   ██ ██    ██ ██    ██ ██      ██   ██ "
echo -e "\t██      ██  ██████  ██████   ██████   ██████  ███████ ██   ██ "                                                                                                                    
echo -e "\t\t\tMulti-Domain Genome Recovery"
echo -e "\t\t\t\tVersion 1.0.0\n\n"
}

help_message () {
        echo""
        echo "Mudoger v=$VERSION"
        echo "Usage: mudoger --module module_name --meta metadata_table.tsv -o output_folder [module_options]"
        echo ""
        echo "  preprocess              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
        echo "          read_qc               explanation (only this submodule)"
        echo "          mem_pred   explanation (only this submodule)"
        echo "          assembly              explanation (only this submodule)"
        
        echo "  prokaryotes              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
        echo "          initial_binning       explanation (only this submodule)"
  
        echo "  viruses              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
               
        echo "  eukaryotes              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
     
        echo "  abundance_tables              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
    
        echo ""
        echo "  --help | -h             show this help message"
        echo "  --version | -v  show metaWRAP version"
        echo "  --show-config   show where the metawrap configuration files are stored"
        echo "";}
  
mudoger_lettering

num_cores=20
megahit=""
metaspades=""

running_location="$(pwd)"

###################################### while
while true; do
	case "$1" in
		--module) active_module=$2; shift 2;;
		--meta) metadata_table=$2; shift 2;;
		-o) output_folder=$2; shift 2;;
		-t) num_cores=$2; shift 2;;
		#-m) memory=$2; shift 2;;
		--metaspades) metaspades="--metaspades"; shift 1;;
		-h | --help) help_message; exit 1; shift 1;;
		--) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done

#############################  


### check if there was any problem with metadata file
out="$(python MuDoGeR/tools/mdcheck.py $metadata_table)"
echo -e "$out"
if [[ $out =~ "Closing" ]]; then
   echo -e "\n TIP: run \"mudoger -h\" for help "
   exit 0
fi
#############################
# if everything is okay, proceeds with asking user for output folder
mkdir -p "$output_folder"


if [ "$active_module" = preprocess ]; then
	echo mudoger preprocess ${@:2}
	module_script=MuDoGeR/src/modules/mudoger-module-1.sh
        ###### loop around samples and run module 1
	aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
	for i in $aux; 
	do 
	r1="$(cat "$metadata_table" | awk -F '\t' '{ if ($1 == "'$i'") {print} }' | cut -f2 | grep '_1.f')"; 
	r2="$(cat "$metadata_table" | awk -F '\t' '{ if ($1 == "'$i'") {print} }' | cut -f2 | grep '_2.f')";
	time $module_script -1 $r1 -2 $r2 -o $output_folder/$i -t $num_cores;
	#echo $module_script -1 $r1 -2 $r2 -o $output_folder/$i -t $num_cores;
	#echo 'test'
	done
	
elif [ "$active_module" = prokaryotes ]; then
	echo mudoger prokaryotes ${@:2}
	module_script=MuDoGeR/src/modules/mudoger-module-2.sh
        ###### loop around samples and run module 1
	aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
	for i in $aux; 
	do 
	r1="$(cat "$metadata_table" | awk -F '\t' '{ if ($1 == "'$i'") {print} }' | cut -f2 | grep '_1.f')"; 
	r2="$(cat "$metadata_table" | awk -F '\t' '{ if ($1 == "'$i'") {print} }' | cut -f2 | grep '_2.f')";
	time  $module_script -1 $r1 -2 $r2 -a $output_folder/$i/assembly/final_assembly.fasta -o $output_folder/$i -t $num_cores;
	done
	
	
	##################################### 
	# 1) GROUP ALL RESULTS FROM ALL LIBRARIES INTO ONE INFORMATIVE TABLE
	# TABLE HAS: BINNAME, LIBRARY, N50, L50, TAXONOMY, COMPLETENESS, CONTAMINATION
	
	# 2) RUN GOTUPICK AND LEAVE RESULTS INSIDE ~output/results/prokaryotes/gotupick
	
	# 3) prepare inputs for BRATS
	
	
	
	
	
	
	
elif [ "$active_module" = viruses ]; then
	echo mudoger viruses ${@:2}
	module_script=MuDoGeR/src/modules/mudoger-module-3.sh
        ###### loop around samples and run module 1
	aux="$(while read l ; do echo "$l" | cut -f1; done < "$metadata_table"  | tr '\t' '\n' | sort |  uniq)";
	for i in $aux; 
	do 
	r1="$(cat "$metadata_table" | awk -F '\t' '{ if ($1 == "'$i'") {print} }' | cut -f2 | grep '_1.f')"; 
	r2="$(cat "$metadata_table" | awk -F '\t' '{ if ($1 == "'$i'") {print} }' | cut -f2 | grep '_2.f')";
	time  $module_script -1 $r1 -2 $r2 -a $output_folder/$i/assembly/final_assembly.fasta -o $output_folder/$i -t $num_cores;
	#echo $module_script -1 $r1 -2 $r2 -a $output_folder/$i/assembly/final_assembly.fasta -o $output_folder/$i -t $num_cores;
	#echo 'test'
	done
else
        comm "Please select a proper module of MuDoGeR."
        help_message
        exit 1
fi


### wraping up prokaryotic results
cd $output_folder
c=1
#### renaming bins and creating map table
mkdir -p results/prokaryotes/otus;
mkdir -p results/prokaryotes/mags;
for bin in  */prokaryotes/binning/unique_bins/*fa;
do nbin=OTU_"$c".fa; 
cp $bin results/prokaryotes/otus/$nbin;
cp $bin results/prokaryotes/mags;
echo -e "$(echo "$bin" | rev | cut -f1 -d'/' | rev | cut -f1 -d'-' )\t\c";  echo $nbin;
c=$((c + 1));
done > results/prokaryotes/map_otus.tsv
### creating metrics information (bbtools)
cat */prokaryotes/metrics/*prok_genomes* | sort | uniq | tac >  results/prokaryotes/genomic_metrics.tsv
### checkm
cat */prokaryotes/metrics/checkm_qc/outputcheckm.tsv | sort | uniq  | grep -v lineage > results/prokaryotes/checkm.tsv
### gtdbtk
cat */prokaryotes/metrics/GTDBtk_taxonomy/*.sum* | sort | uniq  | grep -v 'user_genome' > results/prokaryotes/gtdbtk.tsv

#echo -e "bins\tncontigs" > results/prokaryotes/bins_ncontigs.tsv

echo -e "bins\tncontigs" > results/prokaryotes/header
for mag in results/prokaryotes/mags/*;
do echo -e "$(echo "$mag" | rev | cut -f1 -d'/' | rev)\t\c"; cat $mag | grep ">" | wc -l;
done > results/prokaryotes/aux_bins_ncontigs.tsv  

tac results/prokaryotes/aux_bins_ncontigs.tsv > results/prokaryotes/aux2_bins_ncontigs.tsv
#cat  results/prokaryotes/header results/prokaryotes/aux2_bins_ncontigs.tsv > results/prokaryotes/bins_ncontigs.tsv
paste results/prokaryotes/bins_ncontigs.tsv results/prokaryotes/genomic_metrics.tsv > results/prokaryotes/genomic_metrics_aux.tsv
tac results/prokaryotes/genomic_metrics_aux.tsv > results/prokaryotes/genomic_metrics_aux2.tsv
cat  results/prokaryotes/header results/prokaryotes/genomic_metrics_aux2.tsv > results/prokaryotes/genomic_metrics.tsv



#paste results/prokaryotes/bins_ncontigs.tsv

### run OTU pick
#gOTUpick.sh --bb-input path/to/BBMap-input --checkm-input path/to/CheckM-input --gtdb-input path/to/gtdb-input  -m path/to/mags -o path/to/outputdir

#conda activate gOTUpick
##gOTUpick.sh --bb-input results/prokaryotes/genomic_metrics.tsv --checkm  results/prokaryotes/checkm.tsv --gtdb-input results/prokaryotes/gtdbtk.tsv -m results/prokaryotes/mags -o results/prokaryotes/gOTUpick
#conda deactivate

#conda activate metawrap-env
#echo 'ronaldo'
cd $running_location


script=MuDoGeR/src/scripts/gotupick.sh
chmod +x $script

time $script






###### metawrap configuration
exit 0
#config_file=$(which config-metawrap)
#source $config_file #### DATABASES!!!
#if [[ $? -ne 0 ]]; then
#        echo "cannot find config-metawrap file - something went wrong with the installation!"
#        exit 1
#fi
#if [ "$1" = module_1 ]; then
#        echo mudoger module_1 ${@:2}
#        time ${PIPES}/mudoger-module-1.sh ${@:2}
#        exit 0
#elif [ "$1" = read_qc ]; then
#        time ${PIPES}/read_qc.sh ${@:2}
#        exit 0
#elif [ "$1" = mem_pred ]; then
#        time ${PIPES}/mem_pred.sh ${@:2}
#        exit 0
#elif [ "$1" = assembly ]; then
#        time ${PIPES}/assembly.sh ${@:2}
#        exit 0
#        #### ADD THE REMAINING MODULES AND SUBMODULES
#else
#        comm "Please select a proper module of metaWRAP."
#        help_message
#        exit 1
#fi

