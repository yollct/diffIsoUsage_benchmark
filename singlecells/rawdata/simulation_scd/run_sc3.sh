#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=./unsc3.out
#SBATCH --error=./unsc3.err
#SBATCH --job-name=sc3
#SBATCH --mem-per-cpu=200G
#SBATCH --time=4-00:00:00
#SBATCH --partition=exbio-cpu

set +eu 
eval "$(conda shell.bash hook)"
conda activate /nfs/scratch/chit/.conda/env/r4.3
set -eu
# for i in $(seq 1 3);
# do
#     echo $i
#     R CMD BATCH balance_mat_sim.R
#     mv sim_data_balance/ sim_data_balance${i}_0/
# done

for i in $(seq 1 2);
do
    echo $i
    R CMD BATCH unbalance_mat_sim.R
    mv sim_data_unbalance/ sim_data_unbalance${i}_0/
done

# R CMD BATCH ../metacells/dtu.R
# R CMD BATCH run_dtu_unbalance.R
# R CMD BATCH balance_dropout.R
# R CMD BATCH unbalance_dropout.R
# srun -p exbio-interactive -w norm.exbio.wzw.tum.de --mpi=pmi2 --cpus-per-task=12 --ntasks-per-node=1 --mem=300GB --pty bash