
########## 2 PROKARYOTIC REFINEMENT FOR ARCHAEA (40,30)  ###################
## load conda metawrap conda quality control
conda activate metawrap-env


# arguments declaration
#log="log_qc"                      # definition of path to QC log       

output_folde=$1                # output folder to be created inside master output folder
cores=$2
assembly=$3
con=$4
max=$5
met=$6


metawrap bin_refinement -o "$output_folder" -t $cores -A "$con" -B "$met" -C "$max" -c 40 -x 30

conda deactivate
