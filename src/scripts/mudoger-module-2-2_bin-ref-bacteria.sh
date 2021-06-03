########## 2 PROKARYOTIC REFINEMENT FOR BACTERIA (50,10  ###################
## load conda metawrap conda quality control
conda activate metawrap-env
conda activate metawrap-env

# arguments declaration
log="log_qc"                      # definition of path to QC log       
assembly=$1
forward_library = $2              # forward library path
reverse_library = $3              # reverse library path
output_folder = $4                # output folder to be created inside master output folder
num_cores = 1                     # number of threads


mkdir $prok/bact_ref
outb="$prok/bact_ref"
mkdir $prok/arch_ref
outa="$prok/arch_ref"

con="$prok/concoct_bins"
met="$prok/metabat2_bins"
max="$prok/maxbin2_bins"


metawrap bin_refinement -o "$outb" -t $num_cores -A "$con" -B "$met" -C "$max" -c 50 -x 10
