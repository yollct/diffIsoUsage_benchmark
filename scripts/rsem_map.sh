#!/bin/bash 
# things you don't usually change-------------------

#SBATCH --ntasks=1

#SBATCH --output=/nfs/home/students/chit/output.out

#SBATCH --error=/nfs/home/students/chit/error.err


# things that depend on your task.-----------------

#SBATCH --job-name=rsem


#SBATCH --mem=50G                                      # this is the max capacity of each node in the cluster

#SBATCH --cpus-per-task=2                         # I personally like leaving 2 cpus for node/slurm overheads and then 40 is max

#SBATCH --time=10:00:00
# -------------------------

PATH=$PATH:/nfs/home/students/chit/RSEM/bin/
dir="/nfs/home/students/chit/is_benchmark"

cd ${dir}/alignments

for file in $(ls) 
do
  cd $file
  for sample in $(ls | grep "Transcriptome.out.bam$"); 
  do      
      rsem-calculate-expression --bam --no-bam-output -p 12 --paired-end --forward-prob 0 \
            ${file}Aligned.toTranscriptome.out.bam \
            ${dir}/rsem_ref/asim \
            ${dir}/alignments/${file}/quant 
  done
  cd ..
done
cd ../..