
################### 3 ASSEMBLY  #########################################################
## load conda metawrap conda quality control
conda activate metawrap-env

assembly_folder="$output_master"/ASSEMBLY
mkdir -p "$assembly_folder"      # create assembly output folder

# arguments declaration
log="log_assembly"                # definition of path to assembly log       
forward_library = $1              # forward library path
reverse_library = $2              # reverse library path
output_folder = $3                # output folder to be created inside master output folder
num_cores = 1                     # number of threads
assembler="--metaspades"          # --metaspades or --megahit
# assembly command
metawrap assembly -1 $forward_library -2 $reverse_library -m $mem_gb -t "$num_cores" "$assembler" -o "$assembly_folder"/"$output_folder"

if [ -f "$assembly_folder"/"$output_folder"/final_assembly.fa ]; 
then echo  "assembly was succesful" ; 
else echo "assembly didnt work. please check your resources"; 
fi
