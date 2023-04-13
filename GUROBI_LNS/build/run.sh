#!/bin/bash
#SBATCH --job-name="OP"
#SBATCH --partition=comp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-0:10:00
#SBATCH --array=0-0
##SBATCH --mem-per-cpu=16G

module load cmake/3.15.4-gcc8
module load gurobi/9.0.1

./CSP ${SLURM_ARRAY_TASK_ID}
