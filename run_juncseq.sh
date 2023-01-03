#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=/nfs/scratch/chit/jcseq.out
#SBATCH --error=/nfs/scratch/chit/jcseq.err
#SBATCH --job-name=jcseq-sim
#SBATCH --mem-per-cpu=20G
#SBATCH --time=07-00:00:00

#newgrp docker

#docker run -v /nfs/scratch/lio/isbench_covid/:/MOUNT --rm --name 'jcseq-covid' jcseq
singularity run --bind /nfs/scratch/chit/simulated_real/single_50_8_r1:/MOUNT /nfs/proj/is_benchmark/runner/sing_junctionseq/jcseq.sif 