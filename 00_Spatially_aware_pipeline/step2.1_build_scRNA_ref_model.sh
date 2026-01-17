#!/bin/bash
#SBATCH -J SDAS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16

samp=$(awk -F',' "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" 1.csv)

bash SDAS_build_model.sh $samp
