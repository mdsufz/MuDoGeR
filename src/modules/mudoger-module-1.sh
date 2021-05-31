# help message


forward_library = $1              # forward library path /path/to/libname_1.fastq
reverse_library = $2              # reverse library path. /path/to/libname_2.fastq
output_folder = $3                # master output folder to be created and defined by user


#out_master="$( echo "$output_folder"/"$(echo "$forward_library" | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | cut -f1 -d'_' )")"          # create output master
#mkdir -p $out_master

# 1 QUALITY CONTROL (QC) OF READS
/MuDoGeR/src/scripts/mudoger-module-1-1_QC.sh $forward_library $reverse_library "$3"

# 2 KMER COUNT AND MEMORY ESTIMATION FOR ASSEMBLY
/MuDoGeR/src/scripts/mudoger-module-1-2_kmermempred.sh $forward_library $reverse_library "$3"

# 3 ASSEMBLY
/MuDoGeR/src/scripts/mudoger-module-1-assembly.sh $forward_library $reverse_library "$3"
