#!/bin/bash

# Function to display help
usage() {
    echo "Usage: $0 <module_name> -s <singularity_file_path> -o <output_path> -i <input_data_path> -d <databases_path> -h <home_path> -m <memory> -t <threads> -f <metadata_file> -j <job_name> -N <nodes> -c <cpus> -M <slurm_memory> -T <time_limit> [abundance_tables options]"
    echo "  <module_name>         Module name (e.g., preprocess, prokaryotes, viruses, abundance_tables, eukaryotes)"
    echo "Options:"
    echo "  -s  Path to Singularity file"
    echo "  -o  Path to output folder"
    echo "  -i  Path to input data (metadata table should be here as well)"
    echo "  -d  Path to databases folder"
    echo "  -c  Path for Singularity home directory (for eukaryotes module)"
    echo "  -m  Memory size (for preprocess module)"
    echo "  -t  Number of threads per task for SLURM"
    echo "  -f  Name of the metadata file"
    echo "  -j  Job name for SLURM"
    echo "  -M  Memory per CPU requirement for SLURM (e.g., 1G)"
    echo "  -T  Time limit for SLURM (e.g., 1:00:00 for 1 hour)"
    echo "Abundance Tables Options (only for abundance_tables module):"
    echo "  --reduced (default), --complete, or --genes (exclusive options)"
    echo "  --absolute-values (default), --coverage, and/or --relative-abundance (can be combined)"
    exit 1
}

# Default values
memory=100
thread_count=25
abundance_type="--reduced"  # Default
value_types=("--absolute-values")  # Default as an array
metadata_file="metadata.tsv"  # Default metadata file name
home_path=""

# First argument is the module name
module_name=$1
shift # Shift arguments to parse options

# SLURM parameters defaults
job_name="mudoger_job"
slurm_memory="1G"
time_limit="1:00:00"

# Parsing command line options
while getopts ":s:o:i:d:c:m:t:f:j:M:T:" opt; do
  case $opt in
    s) singularity_file=$(realpath "$OPTARG")
    ;;
    o) output_path=$(realpath "$OPTARG")
    ;;
    i) input_data_path=$(realpath "$OPTARG")
    ;;
    d) databases_path=$(realpath "$OPTARG")
    ;;
    c) home_path=$(realpath "$OPTARG")
    ;;
    m) memory="$OPTARG"
    ;;
    t) thread_count="$OPTARG"
    ;;
    f) metadata_file="$OPTARG"
    ;;
    j) job_name="$OPTARG"
    ;;
    M) slurm_memory="$OPTARG"
    ;;
    T) time_limit="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
       usage
    ;;
  esac
done

# Additional options for abundance_tables
for arg in "$@"; do
  case $arg in
    --reduced|--complete|--genes)
      if [ -n "$exclusive_option_set" ]; then
        echo "Error: --reduced, --complete, and --genes are exclusive options. Only one can be used."
        usage
      else
        abundance_type="$arg"
        exclusive_option_set=true
      fi
      ;;
    --absolute-values|--coverage|--relative-abundance)
      value_types+=("$arg")
      ;;
  esac
done

# Check if all required options are provided
if [ -z "$singularity_file" ] || [ -z "$output_path" ] || [ -z "$input_data_path" ] || [ -z "$databases_path" ]; then
    echo "All options are required."
    usage
fi

# Base command with dynamic paths
base_command="singularity exec --bind $input_data_path:/tools/data_input,$output_path:/tools/mudoger_output_in_container,$databases_path:/tools/dbs"

# Additional Singularity options for eukaryotes module
if [ "$module_name" = "eukaryotes" ] && [ -n "$home_path" ]; then
    base_command="$base_command,$home_path:/mudoger_home --home /mudoger_home"
fi

base_command="$base_command $singularity_file mudoger"

# Construct the command based on the module name
case $module_name in
    preprocess)
        command="$base_command --module preprocess --meta /tools/data_input/$metadata_file -o /tools/mudoger_output_in_container -t $thread_count -m $memory"
        ;;
    prokaryotes|viruses)
        command="$base_command --module $module_name --meta /tools/data_input/$metadata_file -o /tools/mudoger_output_in_container -t $thread_count"
        ;;
    abundance_tables)
        command="$base_command --module abundance_tables --meta /tools/data_input/$metadata_file -o /tools/mudoger_output_in_container -t $thread_count $abundance_type ${value_types[@]}"
        ;;
    eukaryotes)
        if [ -z "$home_path" ]; then
            echo "Home path is required for the eukaryotes module."
            usage
        fi
        command="$base_command --module eukaryotes --meta /tools/data_input/$metadata_file -o /tools/mudoger_output_in_container -t $thread_count"
        ;;
    *)
        echo "Invalid module name"
        usage
        ;;
esac

# SLURM job script creation
echo "Creating slurm Job script!"
job_script="slurm_${job_name}.sh"
echo "#!/bin/bash" > $job_script
echo "#SBATCH --job-name=$job_name" >> $job_script
echo "#SBATCH --cpus-per-task=$thread_count" >> $job_script
echo "#SBATCH --mem-per-cpu=$slurm_memory" >> $job_script
echo "#SBATCH --time=$time_limit" >> $job_script
echo "#SBATCH -o /work/%u/%x-%j.out" >> $job_script
echo "#SBATCH -e /work/%u/%x-%j.err" >> $job_script
echo "SLURM_CPUS_PER_TASK=$thread_count" >> $job_script
echo 'export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}' >> $job_script
echo "$command" >> $job_script

# Submit job script to SLURM
echo "Submitting job to SLURM..."
sbatch $job_script
