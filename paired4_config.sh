#!/bin/bash




##rather to run simulation
sim="false"
##path to the main script
path="/nfs/proj/is_benchmark"
##
#name="pair_50_8_r2_0.5"

meta="/nfs/scratch/chit/new_simulations/${name}/meta.txt"
gtf="/nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.98.gtf"
fasta="/nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
genepred="/nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.98.refflat"
transcript_fasta="/nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.cdna.all.fa"
readfilesdir="/nfs/scratch/chit/new_simulations/${name}/fastq_sim/"
outputdir="/nfs/scratch/chit/new_simulations/${name}/"
index="/nfs/scratch/chit/ref/ens98_star_rsem"


###
depth="40"
nCores="2"
pattern="SRR"
