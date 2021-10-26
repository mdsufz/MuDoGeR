

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
        echo "Usage: mudoger [module]  --meta metadata_table.tsv -o output_folder [module_options]"
        echo ""
        echo "  module_1              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
        echo "          read_qc               explanation (only this submodule)"
        echo "          mem_pred   explanation (only this submodule)"
        echo "          assembly              explanation (only this submodule)"
        
        echo "  module_2              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
        echo "          initial_binning       explanation (only this submodule)"
  
        echo "  module_3              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
               
        echo "  module_4              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
     
        echo "  module_5              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
    
        echo ""
        echo "  --help | -h             show this help message"
        echo "  --version | -v  show metaWRAP version"
        echo "  --show-config   show where the metawrap configuration files are stored"
        echo "";}
  
mudoger_lettering


###################################### while
while true; do
	case "$1" in
		--meta) metadata_table=$2; shift 2;;
		#-2) reverse_library=$2; shift 2;;
		#-o) output_folder=$2; shift 2;;
		#-t) num_cores=$2; shift 2;;
		#-m) memory=$2; shift 2;;
		#--metaspades) metaspades="--metaspades"; shift 1;;
		#-h | --help) help_message; exit 1; shift 1;;
		#--) help_message; exit 1; shift; break ;;
		#*) break;;
	esac
done
############################# end 

echo '--> Please provide a path to the metadata file. ($ mudoger -help -md) for explanation'
read metadata
echo 'path to metadata file: '$metadata
echo 'reading '$metadata
out="$(python MuDoGeR/tools/mdcheck.py $metadata)"
echo -e "$out"
### check if there was any problem with metadata file
if [[ $out =~ "Closing" ]]; then
   echo -e "\n TIP: run \"mudoger -h\" for help "
   exit 0
fi
# if everything is okay, proceeds with asking user for output folder
echo '--> Please provide a path for the output folder. If it does not exist, it will be created'
read output_folder
mkdir -p "$output_folder"
echo '--> Folder created'


md="metadata.tsv" 


if [ "$1" = preprocess ]; then
	echo mudoger preprocess ${@:2}
	module=MuDoGeR/src/modules/mudoger-module-1.sh
        ###### loop around samples and run module 1
	aux="$(while read l ; do echo "$l" | cut -f1; done < "$md"  | tr '\t' '\n' | sort |  uniq)";
	for i in $aux; 
	do 
	r1="$(cat "$md" | awk -F '\t' '{ if ($1 == "'$i'") {print} }' | cut -f2 | grep '_1.f')"; 
	r2="$(cat "$md" | awk -F '\t' '{ if ($1 == "'$i'") {print} }' | cut -f2 | grep '_2.f')";
	time $module -1 $r1 -2 $r2 -o $output_folder/$i ;
	done

else
        comm "Please select a proper module of MuDoGeR."
        help_message
        exit 1
fi


###### metawrap configuration

config_file=$(which config-metawrap)
source $config_file #### DATABASES!!!
if [[ $? -ne 0 ]]; then
        echo "cannot find config-metawrap file - something went wrong with the installation!"
        exit 1
fi
if [ "$1" = module_1 ]; then
        echo mudoger module_1 ${@:2}
        time ${PIPES}/mudoger-module-1.sh ${@:2}
        exit 0
elif [ "$1" = read_qc ]; then
        time ${PIPES}/read_qc.sh ${@:2}
        exit 0
elif [ "$1" = mem_pred ]; then
        time ${PIPES}/mem_pred.sh ${@:2}
        exit 0
elif [ "$1" = assembly ]; then
        time ${PIPES}/assembly.sh ${@:2}
        exit 0
        #### ADD THE REMAINING MODULES AND SUBMODULES
else
        comm "Please select a proper module of metaWRAP."
        help_message
        exit 1
fi

