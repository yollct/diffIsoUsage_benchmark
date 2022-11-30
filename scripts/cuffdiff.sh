#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=mapping.out
#SBATCH --error=mapping.err
#SBATCH --job-name=mapping
#SBATCH --mem-per-cpu=20Gbash 
#SBATCH --time=10:00:00

path="/nfs/home/students/chit/is_benchmark/alignments"
output="/nfs/home/students/chit/is_benchmark/"
PATH=$PATH:/nfs/home/students/chit/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/nfs/home/students/chit/samtools-1.13

cd $path 

# for x in $(ls | grep "sample")
# do
#   cd $x
#   #chmod -R 751 ${x}Aligned.out.bam
#   #samtools sort ${x}Aligned.out.bam > ${x}Aligned.out.sorted.bam
#   samtools view -h -o ${x}Aligned.sortedByCoord.out.sam ${x}Aligned.sortedByCoord.out.bam
#   chmod -R 751 ${x}Aligned.sortedByCoord.out.sam
#   cd ..
# done

group1=$(awk 'BEGIN{FS=OFS=" "}{ if ($2 == 1) print"/nfs/home/students/chit/is_benchmark/alignments/"$1"/"$1"Aligned.sortedByCoord.out.bam" }' ../simulated_reads/sim_rep_info.txt | grep -v "group") 
group2=$(awk 'BEGIN{FS=OFS=" "}{ if ($2 == 2) print"/nfs/home/students/chit/is_benchmark/alignments/"$1"/"$1"Aligned.sortedByCoord.out.bam" }' ../simulated_reads/sim_rep_info.txt | grep -v "group") 
group1l=""
for i in $group1;
do 
    group1l="${group1l},${i}"
done;

group2l=""
for i in $group2;
do 
    group2l="${group2l},${i}"
done;

echo $group1l
cuffdiff -L g1,g2 -u ${outputdir}/merged_asm/merged.gtf -p 20 -b /nfs/data/covid_hscell_4tp/ensembl_99/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa ${group1l:1} ${group2l:1} -o ${output}/results/cuffdiff_results