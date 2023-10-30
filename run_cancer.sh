#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --output=/nfs/scratch/chit/cancer2.out
#SBATCH --error=/nfs/scratch/chit/cancer2.err
#SBATCH --job-name=dtu2
#SBATCH --mem-per-cpu=40G

#path="/nfs/home/students/chit/is_benchmark"

# mkdir ${path}/results
PATH=$PATH:/nfs/home/students/chit/RSEM-1.3.2/


# for dep in 40; do
#     bash main_script.sh -g 30000 -d ${dep} -f 4 -s 2 -r 4 -a /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.cdna.all.fa.gz -o /localscratch/chit/is_sim_fc2
# done;

# bash main_script.sh --config /nfs/data/Sys_CARE/sf_ameling_sys_care/sf_ameling_sys_car540_hgvnvdsx2_podocytes_nephrology/config.sh
bash cancer_main_script.sh --config /nfs/scratch/chit/GSE222260/config.sh