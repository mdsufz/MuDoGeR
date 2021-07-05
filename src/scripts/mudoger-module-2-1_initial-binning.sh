########## 1 INITIAL PROKARYOTIC BINNING  ###################
## load conda metawrap conda quality control
conda activate metawrap-env

# arguments declaration
log="log_qc"                      # definition of path to QC log       
assembly=$1
forward_library = $2              # forward library path
reverse_library = $3              # reverse library path
output_folder = $4                # output folder to be created inside master output folder
num_cores = 1                     # number of threads


#mkdir $prok/bact_ref
#outb="$prok/bact_ref"
#mkdir $prok/arch_ref
#outa="$prok/arch_ref"

#con="$prok/concoct_bins"
#met="$prok/metabat2_bins"
#max="$prok/maxbin2_bins"

metawrap  binning -o "$output_folder" -t "$num_cores" -a "$assembly" --run-checkm --metabat2 --maxbin2 --concoct  "$forward_library" "$reverse_library"
