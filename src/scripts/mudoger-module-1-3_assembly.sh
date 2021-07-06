
################### 3 ASSEMBLY  #########################################################
## load conda metawrap conda quality control
conda activate metawrap-env


# arguments declaration
log="log_assembly"                # definition of path to assembly log       
forward_library=$1                # forward library path
reverse_library=$2                # reverse library path
master_folder=$3                # output folder to be created inside master output folder
num_cores=$4                      # number of threads
assembler="--metaspades"          # --metaspades or --megahit

#out_master="$( echo "$output_folder"/"$(echo "$forward_library" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | cut -f1 -d'_' )")"          # create output master
assembly_folder="$master_folder"/assembly

mkdir -p "$assembly_folder"      # create assembly output folder

kmer_file="$master_folder"/khmer/final_prediction.tsv
mem_mb="$(grep "$(echo "$forward_library" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | cut -f1 -d'_')" $kmer_file | cut -f2)"
 
if [ -z $mem_mb ]; 
then mem_gb=100;
else  mem_gb="$(echo $((mem_mb / 1000)))";
fi
# assembly command
metawrap assembly -1 $forward_library -2 $reverse_library -m $mem_gb -t "$num_cores" "$assembler" -o "$assembly_folder"

if [ -f "$assembly_folder"/final_assembly.fasta ]; 
then echo  "assembly was succesful" ; 
else echo "assembly didnt work. please check your resources"; 
fi
