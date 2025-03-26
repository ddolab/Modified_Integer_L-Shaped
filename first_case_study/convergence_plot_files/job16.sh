#!/bin/bash -l  
#SBATCH --time=6:00:00
#SBATCH --ntasks-per-node=32
#SBATCH --mem=200gb
#SBATCH -p aglarge
#SBATCH --tmp=10gb
module load gurobi
module load julia
julia main_instance_2_16.jl