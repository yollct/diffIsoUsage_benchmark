#!/bin/bash


# things you don't usually change-------------------
#SBATCH --ntasks=1
#SBATCH --output=./star.out
#SBATCH --error=./star.err
# things that depend on your task.-----------------
#SBATCH --job-name=star
#SBATCH --mem=50G                                      # this is the max capacity of each node in the cluster
#SBATCH --cpus-per-task=2                         # I personally like leaving 2 cpus for node/slurm overheads and then 40 is max
#SBATCH --time=20:00:00
# -------------------------
PATH=$PATH:/nfs/home/students/chit/STAR-2.7.9a/bin/Linux_x86_64/
dir="/nfs/home/students/chit/is_benchmark"



for sample in $(ls | grep "fasta$" | awk -F_ '{print $1"_"$2}' | uniq); 
do      

echo $sample
mkdir ./alignments/${sample}

STAR --genomeDir /nfs/data/covid_hscell_4tp/ensembl_106/star_index2.7.9 \
      --runThreadN 30 \
      --readFilesIn  ./simulated_reads/${sample}_1.fasta ./simulated_reads/${sample}_2.fasta \
      --readFilesCommand cat \
      --outSAMstrandField intronMotif \
      --outFilterIntronMotifs RemoveNoncanonical \
      --outSAMattrIHstart 0 \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes jM jI\
      --alignSoftClipAtReferenceEnds No \
      --quantMode GeneCounts TranscriptomeSAM \
      --sjdbGTFfile /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.gtf \
      --sjdbOverhang 100 \
      --outFileNamePrefix ./alignments/${sample}/${sample}
     
done
cd ..

