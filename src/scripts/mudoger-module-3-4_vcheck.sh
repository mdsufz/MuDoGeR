#!/bin/bash


input_uvigs="$1"
output_checkv="$2"
cores="$3"


conda activate checkv-env
export CHECKVDB=/mnt/databases/checkv/checkv-db-v1.0


checkv end_to_end  "$input_uvigs" "$output_checkv" -t "$cores"
