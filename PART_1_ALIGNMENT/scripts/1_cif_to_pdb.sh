#!/bin/bash -l
#SBATCH --cluster="genius"
#SBATCH --account=lp_phylogeo_inf_gpu
#SBATCH --job-name="convert_files"
#SBATCH --ntasks=18
#SBATCH --time=02:00:00
#SBATCH --mem=8G 
#SBATCH --mail-type=END,FAIL           
#SBATCH --mail-user=wout.vanderheyden@student.kuleuven.be

INPUT_DIR="/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/1_cif_files"
OUTPUT_DIR="/data/leuven/358/vsc35887/master_thesis/PCA_af3/TEST_AF_PRED_SIMILARITY/2_pdb_files"

# anders is er een probleem met de storage.
export APPTAINER_CACHEDIR="$VSC_SCRATCH/.apptainer"
export SINGULARITY_CACHEDIR="$VSC_SCRATCH/.apptainer"
export APPTAINER_LIBRARYDIR="$VSC_SCRATCH/.apptainer"
export SINGULARITY_LIBRARYDIR="$VSC_SCRATCH/.apptainer"
export APPTAINERENV_TMPDIR="$VSC_SCRATCH/.apptainer/tmp"
export SINGULARITYENV_TMPDIR="$VSC_SCRATCH/.apptainer/tmp"

## FUNCTION THAT CATCHES THE ERRORS ## 
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

for cif_file in "$INPUT_DIR"/*.cif; do
    base_name=$(basename "$cif_file" .cif)
    pdb_file="$OUTPUT_DIR/$base_name.pdb"

    singularity exec -B "$VSC_DATA:$VSC_DATA" -B "$VSC_SCRATCH:$VSC_SCRATCH" -B "$PWD" \
        docker://kboltonlab/obabel:latest \
        obabel -icif "$cif_file" -opdb -O "$pdb_file"
done

## WRITE CLUSTER METADATA ##
printf END >&2; uptime >&2
trap : 0
echo >&2 '
************
*** DONE ***
************
'
