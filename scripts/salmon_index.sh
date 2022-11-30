#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=mapping.out
#SBATCH --error=mapping.err
#SBATCH --job-name=salmon
#SBATCH --mem=50G
#SBATCH --time=5:00:00
dir="/nfs/home/students/chit/is_benchmark/"

salmon index -t /nfs/data/covid_hscell_4tp/ensembl_106/gentrome.fa.gz -d /nfs/data/covid_hscell_4tp/ensembl_106/decoys.txt -p 30 -i /nfs/data/covid_hscell_4tp/ensembl_106/salmon_index 