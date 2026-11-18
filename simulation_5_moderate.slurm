#!/bin/bash
#SBATCH --job-name=simulation_5_moderate
#SBATCH -o /output
#SBATCH --export=ALL
#SBATCH --partition=cpu-zen4
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL

export R_LIBS="/R_pack"
module load gcc
module load  R/4.4.2-gfbf-2024a

Rscript /simulation_5_moderate.R
