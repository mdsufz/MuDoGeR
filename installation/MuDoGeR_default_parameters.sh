# MuDoGeR default parameters. If the parameter is not explicitly described here, it uses the default paramenters from the original tool. Please read the original tool publications.
#An update is planned to be made to allow editing this file and atomatically adjust MuDoGeR's parameters

# Additional parameters from mudoger-module-1-3_assembly.sh
assembler="--metaspades"  # or "--megahit"

# Additional parameters from mudoger-module-2-2_bin-ref-bacteria.sh
min_completeness="50"
max_contamination="10"

# Parameters from mudoger-module-2-3_bin-ref-archea.sh
min_completeness_archaea="40"
max_contamination_archaea="30"

# Additional parameters from mudoger-module-2-6_bin-QC.sh
tree_type="--reduced_tree"

# Additional parameters from mudoger-module-3-1_viral-investigation.sh
# CLUSTER_GENOMES.pl parameters
coverage_threshold="70"
identity_threshold="95"

# Additional parameters from mudoger-module-3-2_vcontact2.sh
# vcontact2 parameters
cluster_1_path="cluster_one-1.0.jar"

# Additional parameters from mudoger-module-3-3_host-prediction.sh
# WIsH parameters
WIsH_command="predict"
WIsH_bootstrap_mode="1"
WIsH_nullParameters="/path/to/nullParameters.tsv"

# Additional parameters from mudoger-module-4-1_eukrep-eukbin-filter.sh
# Euk bins filter by size
min_euk_size=1500000

# Additional parameters from mudoger-module-4-2_genemark.sh
# Genemark
model="--ES"
min_contig="3000"

# Additional parameters from mudoger-module-4-5_busco.sh
# BUSCO
mode="protein"

# Additional parameters from mudoger-module-5-1_prokOTUpicking.sh
# pOTUpick
fastANI=95

