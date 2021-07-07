



module load Anaconda3/5.3.0
source activate /data/msb/tools/vcheck/conda-vcheck-env
export CHECKVDB=/data/msb/tools/vcheck/checkv-db-v0.6



task=""

date
echo "starting "$task



#command
checkv end_to_end  $1 $2 -t 1
