#!/bin/bash

# things you don't usually change-------------------
#SBATCH --ntasks=1
#SBATCH --output=./htseq.out
#SBATCH --error=./htseq.err

# things that depend on your task.-----------------
#SBATCH --job-name=htseq
#SBATCH --mem=50G                                      # this is the max capacity of each node in the cluster
#SBATCH --cpus-per-task=2                         # I personally like leaving 2 cpus for node/slurm overheads and then 40 is max
#SBATCH --array=1-10
#SBATCH --time=24:00:00
# -------------------------
path="/nfs/home/students/chit/is_benchmark"
PATH=$PATH:/nfs/home/students/chit/samtools-1.13

# mkdir ./results/htseq_exon || true
# for sample in $(ls | grep "sample");
# do  
#     echo $sample
#     #allgroup=$(awk 'BEGIN{FS=OFS=" "}{ print $1"_"$2 }' ${dir}/simulated_reads/sim_rep_info.txt)
#     ngroup=$(awk 'BEGIN{FS=OFS=" "}{ print $1"_"$2 }' ${dir}/simulated_reads/sim_rep_info.txt | grep ${sample} | awk -F_ '{ print $2 }')
#     echo $ngroup
    

#     if [ -f ${dir}/results/htseq_exon/${sample}_${group}_exon.txt ];then
#         echo "${sample} exists"
#         continue
#     else
 
#         samtools view -h ./alignments/${sample}/${sample}Aligned.sortedByCoord.out.bam > ./alignments/${sample}/${sample}Aligned.sortedByCoord.out.sam

#         python ${path}/runner/seqGSEA/count_in_exons.py /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.exon.gtf ./alignments/${sample}/${sample}Aligned.sortedByCoord.out.sam ./results/htseq_exon/${sample}_${group}_exon.txt
#     fi

# done
source /nfs/home/students/chit/is_benchmark/single_config.sh

python ${path}/runner/seqGSEA/count_in_exons.py ${index}/this_gtf.exon ${outputdir}/alignments/SRR12682098/SRR12682098Aligned.sortedByCoord.out.sam ${outputdir} SRR12682098 ${meta}
