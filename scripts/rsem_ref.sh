#!/bin/bash

# things you don't usually change-------------------

#SBATCH --ntasks=1

#SBATCH --output=./rsem.out

#SBATCH --error=./rsem.err

# things that depend on your task.-----------------

#SBATCH --job-name=rsem
#SBATCH --mem=50G                                      # this is the max capacity of each node in the cluster
#SBATCH --array=1-8
#SBATCH --cpus-per-task=2                         # I personally like leaving 2 cpus for node/slurm overheads and then 40 is max

#SBATCH --time=10:00:00
# -------------------------

PATH=$PATH:/nfs/home/students/chit/RSEM/bin/
dir="/nfs/home/students/chit/is_benchmark"

#mkdir ${dir}/rsem_ref
files=("SRR12682097" "SRR12682102" "SRR12682107" "SRR12682112" "SRR12682147" "SRR12682151" "SRR12682155" "SRR12682159")


rsem-prepare-reference --gtf /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.gtf \
/nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${dir}/rsem_ref/asim -p 

for fa in ${files[@]};
do
    rsem-calculate-expression -p 30 \
                            --seed 1234 \
                            --forward-prob 0.5 \
                            --seed-length 25 \
                            --fragment-length-mean -1.0 \
                            --bam ${dir}/simulated_ideal_reads/star/${fa}Aligned.toTranscriptome.out.bam \
                            ${dir}/rsem_ref/asim \
                            ${dir}/simulated_ideal_reads/rsem/sample_${fa} 

    
    # 50 millions reads
    # rsem-simulate-reads ${dir}/rsem_ref/asim \
    #                     ${dir}/simulated_ideal_reads/rsem/sample_${fa}.stat/sample_${fa}.model \
    #                     ${dir}/simulated_ideal_reads/rsem/sample_${fa}.isoforms.results \
    #                     ${dir}/simulated_ideal_reads/rsem/sample_${fa}.stat/sample_${fa}.theta \
    #                     50000000 \
    #                     ${dir}/simulated_ideal_reads/sim_fasta/sample_${fa} 
done

# rsem-calculate-expression --phred64-quals \
#                            -p 8 \
#                            --append-names \
#                            --output-genome-bam \
#                            /nfs/data/covid_hscell_4tp/covid_GSE157490/SRR12682097.unmapped_trimmed.fastq \
#                            ${dir}/rsem_ref/asim \
#                            sample_t0