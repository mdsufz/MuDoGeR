

# help message


forward_library = $1              # forward library path
reverse_library = $2              # reverse library path
output_folder = $3                # output folder to be created inside master output folder





# 1 QUALITY CONTROL (QC) OF READS
/MuDoGeR/src/scripts/mudoger-module-1-1_QC.sh $forward_library $reverse_library $output_folder

# 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY
/MuDoGeR/src/scripts/mudoger-module-1-2_kmermempred.sh $forward_library $reverse_library $output_folder

# 3 ASSEMBLY
