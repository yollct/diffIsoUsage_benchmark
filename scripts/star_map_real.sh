#!/bin/bash


# things you don't usually change-------------------

#SBATCH --ntasks=1

#SBATCH --output=./stareal.out

#SBATCH --error=./stareal.err


# things that depend on your task.-----------------

#SBATCH --job-name=star
#SBATCH --array=1-5

#SBATCH --mem=50G                                      # this is the max capacity of each node in the cluster
#SBATCH --cpus-per-task=2                         # I personally like leaving 2 cpus for node/slurm overheads and then 40 is max
#SBATCH --time=10:00:00
# -------------------------
PATH=$PATH:/nfs/home/students/chit/STAR-2.7.9a/bin/Linux_x86_64/
dir="/nfs/home/students/chit/is_benchmark"
cd ${dir}/simulated_ideal_reads
#uninfected*4, 48hpi*4
files=("SRR12682097" "SRR12682102" "SRR12682107" "SRR12682112" "SRR12682147" "SRR12682151" "SRR12682155" "SRR12682159")
#"SRR12682097" "SRR12682102" "SRR12682107"


for sample in ${files[@]};
do    
      echo $sample
      STAR --genomeDir /nfs/data/covid_hscell_4tp/ensembl_106/star_index2.7.9/ \
            --runThreadN 30 \
            --readFilesIn  /nfs/data/covid_hscell_4tp/covid_GSE157490/${sample}.unmapped_trimmed.fastq \
            --readFilesCommand cat \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMattrIHstart 0 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM MD\
            --alignSoftClipAtReferenceEnds No \
            --quantMode GeneCounts TranscriptomeSAM \
            --sjdbGTFfile /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.gtf \
            --sjdbOverhang 100 \
            --outFileNamePrefix star/${sample} 
done


