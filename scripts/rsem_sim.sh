#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=/nfs/scratch/chit/rsem_sim.out
#SBATCH --error=/nfs/scratch/chit/rsem_sim.err
#SBATCH --job-name=rsem_sim
#SBATCH --mem-per-cpu=60G
#SBATCH --time=10:00:00


PATH=$PATH:/nfs/home/students/chit/RSEM-1.3.3
dir="/nfs/home/students/chit/is_benchmark"

files=("SRR12682098" "SRR12682103" "SRR12682108" "SRR12682113" "SRR12682147" "SRR12682151" "SRR12682155" "SRR12682159") ### SE covid ts
#files=$(cat /nfs/scratch/chit/isbench_covid/simlist.txt)
readfilesdir="/nfs/scratch/chit/covid_ts/"
index="/nfs/scratch/chit/ensembl_106/star_rsem"
gtf="/nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.gtf"
fasta="/nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
outputdir="/nfs/scratch/chit/simulated_real"

! test -f ${outputdir}/fastq_sim || mkdir ${outputdir}/fastq_sim

for fa in ${files[@]}; do
    ## 50 million
    echo ${fa}
    rsem-simulate-reads ${index}/asim \
                            ${readfilesdir}/rsem_out/${fa}/${fa}.stat/${fa}.model \
                            ${outputdir}/rsem_sim/${fa}.isoforms.results \
                            ${readfilesdir}/rsem_out/${fa}/${fa}.stat/${fa}.theta \
                            0.2 \
                            50000000 \
                            ${outputdir}/fastq_sim/${fa}

done