########## 2 PROKARYOTIC REFINEMENT FOR BACTERIA (50,10  ###################
## load conda metawrap conda quality control
conda activate metawrap-env



# arguments declaration
log="log_qc"                      # definition of path to QC log       

output_folder = $1                # output folder to be created inside master output folder
cores=$2
assembly=$3
con=$4
met=$5
max=$6


metawrap bin_refinement -o "$output_folder" -t $cores -A "$con" -B "$met" -C "$max" -c 50 -x 10

conda deactivate

