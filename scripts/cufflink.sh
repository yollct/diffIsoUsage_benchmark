#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=mapping.out
#SBATCH --error=mapping.err
#SBATCH --job-name=cufflink
#SBATCH --array=1-10
#SBATCH --mem-per-cpu=20G
#SBATCH --time=100:00:00
##### this cufflinks uses python 2


path="/nfs/home/students/chit/is_benchmark/results"
output="/nfs/home/students/chit/is_benchmark/"
PATH=$PATH:/nfs/home/students/chit/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/nfs/home/students/chit/samtools-1.13

cd ${output}/alignments

    #samtools sort -O bam -o ${file}.sorted ${file}
for sample in $(ls | grep "sample");
do  
    cd $sample
    
    if [ -f ${path}/cufflinks/${sample}/isoform.fpkm_tracking ]; then
        echo ${sample} exists
        continue
    else
        echo $sample
        mkdir ${path}/cufflinks/${sample}
        cufflinks -g /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.gtf -o ${path}/cufflinks/${sample} ${sample}Aligned.sortedByCoord.out.bam -p 40
    fi
    cd ..
done

cd ${output}/results
find ${output}/results -type f | grep "transcripts.gtf" > assembly_list.txt
/nfs/home/students/chit/cufflinks-2.2.1.Linux_x86_64/cuffmerge assembly_list.txt -ref-gtf /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.chr_patch_hapl_scaff.gtf -s /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.dna.primary_assembly.fa -o ${output}/results



