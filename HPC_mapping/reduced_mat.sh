#!/bin/bash

#SBATCH -q batch
#SBATCH --job-name=temptest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=23:59:59
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c-oconnc@jax.org

## R singularity image file
sif_path="/pod/2/reinholdt-lab/Callan/invitroas/scripts/r-packages-qtl2_1.1.sif"
## R script to distribute
rscript_path="/pod/2/reinholdt-lab/Callan/spring_2023/scripts/scanmat_permsonly.R"
module load singularity
singularity exec $sif_path Rscript $rscript_path
