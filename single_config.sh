#!/bin/bash

##rather to run simulation
sim="false"
##path to the main script
path="/nfs/home/students/chit/is_benchmark"
##
name="single_50_8_r2_0.5"

meta="/nfs/scratch/chit/simulated_real/${name}/meta.txt"
gtf="/nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.98.gtf"
gff="/nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.98.gff3"
fasta="/nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
genepred="/nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.98.refflat"
transcript_fasta="/nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.cdna.all.fa"
readfilesdir="/nfs/scratch/chit/simulated_real/${name}/fastq_sim/"
outputdir="/nfs/scratch/chit/simulated_real/${name}/"
#index="/nfs/scratch/chit/ensembl_106/star_rsem"
index="/nfs/scratch/chit/ref/ens98_star_rsem"
###


###
depth="40"
nCores="2"
pattern="SRR"
