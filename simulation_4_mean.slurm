#!/bin/bash
#SBATCH --job-name=mean_400
#SBATCH -o /output/out
#SBATCH --export=ALL
#SBATCH --partition=cpu-zen4
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem=30G
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL

export R_LIBS="/R_pack"
module load gcc
module load  R/4.4.2-gfbf-2024a

Rscript /simulation_4_mean.R
