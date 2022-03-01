#!/bin/bash

# loading conda environment
echo '------- START MODULE 2-6 BIN Quality - CheckM'
# loading conda environment
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database


conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/metawrap_env



extension="fa"

#Set it only once during database configuration
#checkm data setRoot "$CHECKM_DB"

# necessary path variable
checkm lineage_wf -t $2 --reduced_tree --tab_table -x $extension -f "$1"/metrics/checkm_qc/outputcheckm.tsv $input_bins_folder $output_results


conda deactivate
