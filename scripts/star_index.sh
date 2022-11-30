#!/bin/bash


# things you don't usually change-------------------

#SBATCH --ntasks=1

#SBATCH --output=./staroutput.out

#SBATCH --error=./starerror.err


# things that depend on your task.-----------------

#SBATCH --job-name=star


#SBATCH --mem=50G                                      # this is the max capacity of each node in the cluster

#SBATCH --cpus-per-task=10                         # I personally like leaving 2 cpus for node/slurm overheads and then 40 is max

#SBATCH --time=10:00:00
# -------------------------
PATH=$PATH:/nfs/home/students/chit/STAR-2.7.9a/bin/Linux_x86_64/
dir="/nfs/data/covid_hscell_4tp/ensembl_106"
## for ensembl gtf
# STAR --runMode genomeGenerate --genomeDir star_index \
#             --genomeFastaFiles annotations/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
#             --sjdbGTFfile annotations/Homo_sapiens.GRCh38.99.gtf \
#             --sjdbOverhang 50 --outFileNamePrefix v99
            
## for simulation
rm -rf ${dir}/star_index2.7.9
mkdir ${dir}/star_index2.7.9
STAR --runMode genomeGenerate --genomeDir ${dir}/star_index2.7.9 \
            --genomeFastaFiles /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
            --sjdbGTFfile /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.gtf \
            --sjdbOverhang 100 --outFileNamePrefix v106 --runThreadN 30