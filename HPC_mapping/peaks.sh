#!/bin/bash
#SBATCH --job-name=1_median
#SBATCH -q batch
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --mem 100GB # memory pool for all cores
#SBATCH -t 2-12:00 # time (D-HH:MM)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c-oconnc@jax.org

## R singularity image file
sif_path="/pod/2/reinholdt-lab/Callan/qtl2_175/singledose_correction/r-packages-qtl2_1.1.sif"
## R script to distribute
rscript_path="/pod/2/reinholdt-lab/Callan/spring_2023/scripts/getpeaks.R"

module load singularity
singularity exec $sif_path Rscript $rscript_path

