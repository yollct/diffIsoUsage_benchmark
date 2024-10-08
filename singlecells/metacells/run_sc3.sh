#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=./unsc3.out
#SBATCH --error=./unsc3.err
#SBATCH --job-name=sc3
#SBATCH --mem-per-cpu=200G
#SBATCH --time=4-00:00:00
#SBATCH --partition=exbio-cpu


eval "$(conda shell.bash hook)"

conda activate r-env
echo "conda activated"

# for i in $(seq 1 3);
# do
#     echo $i
#     R CMD BATCH balance_mat_sim.R
#     mv sim_data_balance sim_data_balance${i}_0.1
# done

R CMD BATCH unbalance_dtu.R
# R CMD BATCH run_dtu_unbalance.R
# R CMD BATCH balance_dropout.R
# R CMD BATCH unbalance_dropout.R
# srun -p exbio-interactive -w norm.exbio.wzw.tum.de --mpi=pmi2 --cpus-per-task=12 --ntasks-per-node=1 --mem=300GB --pty bash