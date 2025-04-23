#!/bin/bash
#SBATCH --clusters=genius
#SBATCH --account="lp_phylogeo_inf_gpu"
#SBATCH --partition=dedicated_rega_gpu
#SBATCH --gres=gpu:1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=18
#SBATCH --mem 170000M
#SBATCH --time="168:00:00"

####################################
##### Run a beast analysis ######
####################################

#### How to run: ####
# sbatch beast.slurm <my-xml-file>

# write a function that catches any errors
abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

set -e
trap 'abort' 0

# write some cluster metadata to the error and output files
hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

# Check if an XML file was provided
if [ $# -eq 0 ]; then
    echo "Error: No XML file specified" >&2
    echo "Usage: sbatch beast.slurm <my-xml-file>" >&2
    exit 1
fi

module load cluster/genius/dedicated_rega_gpu
module load CUDA

file=$(basename $1)

echo "Running BEAST on $1" >&2

singularity exec --nv -B $VSC_DATA:$VSC_DATA -B $VSC_STAGING:$VSC_STAGING -B $VSC_SCRATCH:$VSC_SCRATCH -B $VSC_HOME:$VSC_HOME -B $PWD docker://jklaps/beast-beagle-cuda:v10.5.0-beta5 \
    beast \
    -beagle_GPU \
    -beagle_order 1 \
    -beagle_double \
    -working \
    -save_stem $file \
    $1 \
    > $file.beast.out

# write some cluster metadata to the error and output files
printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE ***
************
'