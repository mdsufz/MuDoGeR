

help_message () {
        echo""
        echo "Mudoger v=$VERSION"
        echo "Usage: Mudoger [module]"
        echo ""
        echo "  module-1              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
        echo "  read_qc               explanation (only this submodule)"
        echo "  kmer mem prediction   explanation (only this submodule)"
        echo "  assembly              explanation (only this submodule)"
        
        echo "  module-2              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
        echo "  initial_binning       explanation (only this submodule)"
  
        echo "  module-3              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
               
        echo "  module-4              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
     
        echo "  module-5              runs all steps from module 1 (read_qc, kmer mem prediction and assembly)"
    
     
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
