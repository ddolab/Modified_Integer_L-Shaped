#!/bin/bash -l  
#SBATCH --time=6:00:00
#SBATCH --ntasks-per-node=64
#SBATCH --mem=150gb
#SBATCH --tmp=10gb
module load gurobi
module load julia/1.8.0
julia main_instance.jl $SLURM_ARRAY_TASK_ID