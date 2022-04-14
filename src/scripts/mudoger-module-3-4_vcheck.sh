#!/bin/bash

conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source $config_path
source $database


input_uvigs="$1"
output_checkv="$2"
cores="$3"


conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH"/checkv_env


export CHECKVDB="$DATABASES_LOCATION"checkv/checkv-db-v1.0


checkv end_to_end  "$input_uvigs" "$output_checkv" -t "$cores"
