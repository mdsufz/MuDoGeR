#!/bin/bash
# PROKARYOTES MODULES FILESYSTEM STRUCTURE
#   BINNING
#     INITIAL-BINNING
#     REFINEMENT-BACTERIA
#     REFINEMENT-ARCHAEA
#     BIN-FINAL-SET (REDUNDANCY REMOVAL)
#   METRICS
#     CHECKM
#     GTDBtk
#     PROKKA
#     STATS (GENOMIC)
#       N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more...



libname=$1
assembly = $2
forward_library = $3              # forward library path /path/to/libname_1.fastq
reverse_library = $4             # reverse library path. /path/to/libname_2.fastq
output_folder = $5              # master output folder to be created and defined by user
cores=1


#lib="$( echo "$output_folder"/"$(echo "$forward_library" | rev | cut -f2 -d'/' | rev | cut -f1 -d'.' | cut -f1 -d'_' )")"          # create output master

mkdir -p "$libname"/prokaryotes

mkdir -p "$libname"/prokaryotes/binning

# 1 INITIAL BINNING USING METAWRAP (CONCOCT, METABAT2, MAXBIN2)
MuDoGeR/src/scripts/mudoger-module-2-1_initial-binning.sh "$assembly"          \
                                      "$forward_library"   \
                                      "$reverse_library"   \
                                      "$libname"/prokaryotes/binning/initial-binning


# 2 BIN REFINEMENT USING METAWRAP FOR BACTERIA (50,10)
MuDoGeR/src/scripts/mudoger-module-2-2_bin-ref-bacteria.sh "$libname"/prokaryotes/binning/refinement-bac                \
                                       "$cores"                                                     \
                                       "$assembly"                                                  \
                                       "$libname"/prokaryotes/binning/initial-binning/concoct_bins  \
                                       "$libname"/prokaryotes/binning/initial-binning/maxbin2_bins  \
                                       "$libname"/prokaryotes/binning/initial-binning/metabat2_bins \
              
              

# 3 BIN REFINEMENT USING METAWRAP FOR ARCHAEA (40,30)
MuDoGeR/src/scripts/mudoger-module-2-3_bin-ref-archea.sh "$libname"/prokaryotes/binning/refinement-arc                \
                                     "$cores"                                                     \
                                     "$assembly"                                                  \
                                     "$libname"/prokaryotes/binning/initial-binning/concoct_bins  \
                                     "$libname"/prokaryotes/binning/initial-binning/maxbin2_bins  \
                                     "$libname"/prokaryotes/binning/initial-binning/metabat2_bins \
              
# 4 BIN REDUNANCY REMOVAL
cd "$libname"/prokaryotes/binning

lib=$libname;
bin_count=0; 
mkdir -p unique_bins; 
md5sum refinement-bac/metawrap*bins/*fa  >  unique_bins/md5_sum ;
md5sum refinement-arc/metawrap*bins/*fa  >> unique_bins/md5_sum ;
cat unique_bins/md5_sum | cut -f1 -d' ' | sort | uniq > unique_bins/md5_unique; 

while read l; 
do 
bininit="$(grep "$l" unique_bins/md5_sum  | head -n1 | cut -f3 -d' ')";
binafter=unique_bins/"$lib"-bin."$bin_count".fa; 
cp $bininit $binafter ;
bin_count=$[$bin_count+1];
done < unique_bins/md5_unique ; 


# 5 GTDBtk taxonomy assignment




# 6 CheckM quality control

# 7 PROKKA Annotation

# 8 Metrics: N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more...

# 9 WRAPUP INTO NICE TABLE

#
