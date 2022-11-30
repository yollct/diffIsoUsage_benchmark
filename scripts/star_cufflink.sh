#!/bin/bash


# things you don't usually change-------------------

#SBATCH --ntasks=1

#SBATCH --output=/nfs/home/students/chit/output.out

#SBATCH --error=/nfs/home/students/chit/error.err


# things that depend on your task.-----------------

#SBATCH --job-name=star


#SBATCH --mem=50G                                      # this is the max capacity of each node in the cluster

#SBATCH --cpus-per-task=2                         # I personally like leaving 2 cpus for node/slurm overheads and then 40 is max

#SBATCH --time=10:00:00
# -------------------------

bash star_map.sh
bash cufflink.sh