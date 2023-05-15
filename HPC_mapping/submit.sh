#!/bin/bash

#SBATCH --qos=batch
#SBATCH --partition=compute
#SBATCH --job-name=qtl_mapping
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12gb
#SBATCH --time=2-23:59:59
#SBATCH --array=1-5375%100
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c-oconnc@jax.org

## Command line arguments
# first argument:    viewer_path         -- path to raw protein/peptide data
# second argument:   dataset             -- dataset list object in QTL Viewer
# third argument:    data_type           -- Form of data to be used: norm or rz
# fourth argument:   output_dir          -- path to write output directory
# fifth argument:    use_kinship         -- logical on whether to use kinship
# sixth argument:    num                 -- number of NAs or 0s
# seventh argument:  key_value           -- NA, zero, or duplicate
# eighth argument:   this_chunk          -- chunk out of 100

## R singularity image file
sif_path="/pod/2/reinholdt-lab/Callan/invitroas/scripts/r-packages-qtl2_1.1.sif"
## R script to distribute
rscript_path="/pod/2/reinholdt-lab/Callan/spring_2023/scripts/qtlmapping.R"
module load singularity
singularity exec $sif_path Rscript $rscript_path viewer_path=/pod/2/reinholdt-lab/Callan/perms2/scripts/qtldata.Rdata dataset=/pod/2/reinholdt-lab/Callan/spring_2023/input/pheno_rz.rds output_dir=/pod/2/reinholdt-lab/Callan/spring_2023/output/ use_kinship=TRUE num=40 key_value=NA this_chunk=$SLURM_ARRAY_TASK_ID
