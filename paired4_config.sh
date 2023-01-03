#!/bin/bash

##rather to run simulation
sim="false"
##path to the main script
path="/nfs/home/students/chit/is_benchmark"
##
name="pair_50_8_r2"

meta="/nfs/scratch/chit/simulated_real/${name}/meta.txt"
gtf="/nfs/data/references/ensembl98_GRCh38/Homo_sapiens.GRCh38.98.gtf"
fasta="/nfs/data/references/ensembl98_GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
transcript_fasta="/nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.cdna.all.fa"
readfilesdir="/nfs/scratch/chit/simulated_real/${name}/fastq_sim/"
outputdir="/nfs/scratch/chit/simulated_real/${name}/"
index="/nfs/scratch/chit/ref/ens98_star_rsem"


###
depth="40"
nCores="16"
pattern="SRR"
