# PROKARYOTES
# BE AWARE OF QUOTE MARKS "". THEY MIGHT BECOME PROBLEMATIC


mkdir -p "$output_master"       # create master output folder


if [[ -d "$3/RAW_READS/$l" ]]; then
echo "$l already in progress/exists"
continue
else 


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


mkdir $prok/bact_ref
outb="$prok/bact_ref"
mkdir $prok/arch_ref
outa="$prok/arch_ref"

con="$prok/concoct_bins"
met="$prok/metabat2_bins"
max="$prok/maxbin2_bins"

metawrap  binning -o $bin -t ${NSLOTS:-1} -a $assembly --run-checkm --metabat2 --maxbin2 --concoct  $forward_read $reverse_read

metawrap bin_refinement -o "$outb" -t $num_cores -A "$con" -B "$met" -C "$max" -c 50 -x 10

metawrap bin_refinement -o "$outa" -t $num_cores -A "$con" -B "$met" -C "$max" -c 40 -x 30

ch_inpb="$outb/metawrap_50_10_bins"
ch_inpa="$outa/metawrap_40_30_bins"


metawrap_binning "$2" "$3" "$4" "$5" )

/metawrap_bin_refinement "$2"/refinement-50-10 "$2"/concoct_bins "$2"/maxbin2_bins "$2"/metabat2_bins

/metawrap_bin_refinement-40-30.sh "$2"/refinement-40-30 "$2"/concoct_bins "$2"/maxbin2_bins "$2"/metabat2_bins



