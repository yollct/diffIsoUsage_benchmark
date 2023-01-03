#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=/nfs/scratch/chit/paired8.out
#SBATCH --error=/nfs/scratch/chit/paired8.err
#SBATCH --job-name=paired8
#SBATCH --mem-per-cpu=20G
#SBATCH --time=07-00:00:00

#path="/nfs/home/students/chit/is_benchmark"

# mkdir ${path}/results


# for dep in 40; do
#     bash main_script.sh -g 30000 -d ${dep} -f 4 -s 2 -r 4 -a /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.cdna.all.fa.gz -o /localscratch/chit/is_sim_fc2
# done;

# bash main_script.sh --config /nfs/data/Sys_CARE/sf_ameling_sys_care/sf_ameling_sys_car540_hgvnvdsx2_podocytes_nephrology/config.sh
#bash main_script.sh --config /nfs/home/students/chit/is_benchmark/paired_config.sh
bash main_script.sh --config /nfs/home/students/chit/is_benchmark/paired4_config.sh
# source /nfs/home/students/chit/is_benchmark/paired_config.sh
# if ! test -f ${index}/salmon_index;
#   then mkdir -p ${index}/salmon_index || true
#     echo "Building Salmon index"
#     grep '^>' <(cat ${fasta}) | cut -d ' ' -f 1 > ${index}/decoys.txt
#     sed -i.bak -e 's/>//g' ${index}/decoys.txt
#     cat ${transcript_fasta} ${fasta} > ${index}/gentrome.fa.gz
#     salmon index -t ${index}/gentrome.fa.gz -d ${index}/decoys.txt -p 12 -i ${index}/salmon_index 
#   else echo "Reusing Salmon index"
# fi 
# PATH=$PATH:/nfs/home/students/chit/tophat-2.1.1.Linux_x86_64
# PATH=$PATH:/nfs/home/students/chit/bowtie2-2.4.5-linux-x86_64
#source /nfs/home/students/chit/is_benchmark/single_config.sh

# tophat2 --GTF ${gtf} --no-novel-juncs --no-novel-indels \
# 	  -r 800 -p 10 --no-coverage-search -o ${outputdir}/alignments/tophat/{}/ --transcriptome-index ${bowtie_index}/transcript_index/trans_in \
# 	  ${bowtie_index}/genome_index/genome_in \
#     ${readfilesdir}/sample_01_1.fasta.gz ${readfilesdir}/sample_01_2.fasta.gz

# rsem-calculate-expression \
#             --paired-end \
#             --star \
#             --star-path /nfs/home/students/chit/STAR-2.7.4a/bin/Linux_x86_64 \
#             --star-gzipped-read-file \
#             -p 12 \
#             /nfs/proj/is_benchmark/simulated_reads/sample_01_1.fasta.gz \
#             /nfs/proj/is_benchmark/simulated_reads/sample_01_2.fasta.gz \
#             /nfs/data/covid_hscell_4tp/ensembl_106/star_index2.7.9/asim \
#             /nfs/proj/is_benchmark/rsem_out/sample_01/sample_01

