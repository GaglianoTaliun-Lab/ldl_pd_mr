#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm_manhattan.log

module load r/4.0.5
Rscript run_manhattan.r
