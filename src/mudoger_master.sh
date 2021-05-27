



mudoger_header () {
        echo "                #~  /$$      /$$           /$$$$$$$             /$$$$$$            /$$$$$$$ "
        echo "                #~ | $$$    /$$$          | $$__  $$           /$$__  $$          | $$__  $$"
        echo "                #~ | $$$$  /$$$$ /$$   /$$| $$  \ $$  /$$$$$$ | $$  \__/  /$$$$$$ | $$  \ $$"
        echo "                #~ | $$ $$/$$ $$| $$  | $$| $$  | $$ /$$__  $$| $$ /$$$$ /$$__  $$| $$$$$$$/"
        echo "                #~ | $$  $$$| $$| $$  | $$| $$  | $$| $$  \ $$| $$|_  $$| $$$$$$$$| $$__  $$"
        echo "                #~ | $$\  $ | $$| $$  | $$| $$  | $$| $$  | $$| $$  \ $$| $$_____/| $$  \ $$"
        echo "                #~ | $$ \/  | $$|  $$$$$$/| $$$$$$$/|  $$$$$$/|  $$$$$$/|  $$$$$$$| $$  | $$"
        echo "                #~ |__/     |__/ \______/ |_______/  \______/  \______/  \_______/|__/  |__/" }																								 



help_message () {
        echo""
        echo "Mudoger v=$VERSION"
        echo "Usage: Mudoger [module]"
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







