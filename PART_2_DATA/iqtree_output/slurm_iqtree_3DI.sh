#!/bin/bash
#SBATCH --account="lp_phylogeo_inf_gpu"
#SBATCH --job-name="master_thesis"
#SBATCH --cluster="genius"
#SBATCH --partition="batch"
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH -t 01:30:00
#SBATCH --mem=5G
#SBATCH --mail-type=END,FAIL           
#SBATCH --mail-user=wout.vanderheyden@student.kuleuven.be
#SBATCH --chdir=/data/leuven/358/vsc35887/master_thesis/data_flaviviridae/glycoprotein_structural_alignments_and_trees/3di

export VSC_DATA="/data/leuven/358/vsc35887"

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

## RUN 3DI ##
singularity exec --nv -B $VSC_DATA:$VSC_DATA  -B $PWD docker://quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0 \
    iqtree2 -s refolded_fullglyco_E2_3di_famsa.fas \
    -m 3DI -bb 1000 \
    -alrt 1000 \
    -nt AUTO \
    -pre $VSC_DATA/master_thesis/iqtree_output/refolded_fullglyco_E2_famsa_iqtree/refolded_fullglyco_E2_3di_famsa

## RUN 3DI trim35##
singularity exec --nv -B $VSC_DATA:$VSC_DATA  -B $PWD docker://quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0 \
    iqtree2 -s refolded_fullglyco_E2_3di_famsa_trim35.fas \
    -m 3DI -bb 1000 \
    -alrt 1000 \
    -nt AUTO \
    -pre $VSC_DATA/master_thesis/iqtree_output/refolded_fullglyco_E2_famsa_iqtree/refolded_fullglyco_E2_3di_famsa_trim35

## RUN AA ##
singularity exec --nv -B $VSC_DATA:$VSC_DATA -B $PWD docker://quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0 \
    iqtree2 -s refolded_fullglyco_E2_AA_famsa.fas \
    -m TEST -bb 1000 \
    -nt AUTO \
    -pre $VSC_DATA/master_thesis/iqtree_output/refolded_fullglyco_E2_famsa_iqtree/refolded_fullglyco_E2_AA_famsa

## RUN AA trim35##
singularity exec --nv -B $VSC_DATA:$VSC_DATA -B $PWD docker://quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0 \
    iqtree2 -s refolded_fullglyco_E2_AA_famsa_3ditrim35.fas \
    -m TEST -bb 1000 \
    -nt AUTO \
    -pre $VSC_DATA/master_thesis/iqtree_output/refolded_fullglyco_E2_famsa_iqtree/refolded_fullglyco_E2_AA_famsa_3ditrim35


# write some cluster metadata to the error and output files
printf END >&2; uptime >&2
trap : 0
echo >&2 '
************
*** DONE ***
************
'