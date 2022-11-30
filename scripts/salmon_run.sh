#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=mapping.out
#SBATCH --error=mapping.err
#SBATCH --job-name=mapping
#SBATCH --mem-per-cpu=20G
#SBATCH --time=24:00:00

path="/nfs/home/students/chit/is_benchmark/simulated_reads"
output="/nfs/home/students/chit/is_benchmark/"
PATH=$PATH:/nfs/home/students/chit/.conda/envs/salmon/bin

#uninfected*4, 48hpi*4
#files=("SRR12682102" "SRR12682107" "SRR12682112" "SRR12682147" "SRR12682151" "SRR12682155" "SRR12682159")
#"SRR12682097"




for file in $(ls | grep "fasta$" | awk -F_ '{print $1"_"$2}' | uniq); 
do
        #nsam=$(ls | grep $file | grep "fastq.gz" | grep -v "md5\|zip" | wc -l)
        mkdir ./salmon_out/${file} || true
        if [ -f ./salmon_out/${file}/quant.sf ]; then
                echo ${file} exists 
                continue
        else
                sample=$file
                echo $file

                salmon quant -p 20 -i /nfs/data/covid_hscell_4tp/ensembl_106/salmon_index -l IU -1 ./simulated_reads/${file}_1.fasta  -2 ./simulated_reads/${file}_2.fasta --validateMappings -o ./salmon_out/${file} 
                #salmon quant -t ${output}/annotations/splicing_variants.fa -l IU -a ${file}Aligned.toTranscriptome.out.bam -o ${output}/salmon_out/${file} -p 5
        fi
done
