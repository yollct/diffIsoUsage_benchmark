#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=/nfs/scratch/chit/drim81.out
#SBATCH --error=/nfs/scratch/chit/drim81.err
#SBATCH --job-name=drim81
#SBATCH --mem-per-cpu=30G
#SBATCH --time=4-00:00:00
#SBATCH --partition=exbio-cpu

#path="/nfs/home/students/chit/is_benchmark"

# mkdir ${path}/results


# for dep in 40; do
#     bash main_script.sh -g 30000 -d ${dep} -f 4 -s 2 -r 4 -a /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.cdna.all.fa.gz -o /localscratch/chit/is_sim_fc2
# done;


seq="pair"
nsam="8"
bg="0.1"
reps="r4 r5"

for rep in $reps;
    do
        Rscript /nfs/proj/is_benchmark/runner/rerun_drimseq.R /nfs/scratch/chit/new_simulations/${seq}_50_${nsam}_${rep}_${bg} /nfs/proj/is_benchmark /nfs/scratch/chit/new_simulations/${seq}_50_${nsam}_${rep}_${bg}/meta.txt "N" "T"
    done;
#         done;

