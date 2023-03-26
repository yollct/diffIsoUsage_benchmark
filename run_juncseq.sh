#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=/nfs/proj/is_benchmark/jcseq.out
#SBATCH --error=/nfs/proj/is_benchmark/jcseq.err
#SBATCH --job-name=simjcseq
#SBATCH --mem-per-cpu=20G
#SBATCH --time=07-00:00:00

#newgrp docker

#docker run -v /nfs/scratch/chit/simulated_real/single_50_8_r1:/MOUNT --rm --name 'jcseq-covid' jcseq
singularity run --bind /nfs/scratch/chit/simulated_real/single_50_8_r1:/MOUNT --bind /nfs/proj/is_benchmark/Rlib:/usr/local/lib/R/site-library /nfs/proj/is_benchmark/runner/sing_junctionseq/jcseq.sif 
