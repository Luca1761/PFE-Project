#!/bin/bash
#SBATCH --array=0-149
#SBATCH --cpus-per-task=16
#SBATCH --mem=36G
#SBATCH --time=06:30:00
#SBATCH --job-name=expes_DSIRP
#SBATCH -o expes_out/expes_%A.out
#SBATCH -e expes_out/expes_%A.errarray

python3.10 expes.py --expe_id=$SLURM_ARRAY_TASK_ID