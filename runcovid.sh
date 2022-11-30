#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=/nfs/scratch/chit/covid_d.out
#SBATCH --error=/nfs/scratch/chit/covid_d.err
#SBATCH --job-name=covid_d
#SBATCH --mem-per-cpu=350 G
#SBATCH --time=07-00:00:00


#path="/nfs/home/students/chit/is_benchmark"
PATH=$PATH:/nfs/home/students/chit/.conda/envs/nease/bin
# mkdir ${path}/results


# for dep in 40; do
#     bash main_script.sh -g 30000 -d ${dep} -f 4 -s 2 -r 4 -a /nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.cdna.all.fa.gz -o /localscratch/chit/is_sim_fc2
# done;

# bash main_script.sh --config /nfs/data/Sys_CARE/sf_ameling_sys_care/sf_ameling_sys_car540_hgvnvdsx2_podocytes_nephrology/config.sh
bash main_script.sh --config /nfs/home/students/chit/is_benchmark/covid_config.sh
#source /nfs/home/students/chit/is_benchmark/covid_config.sh
#bash downsample_script.sh --config /nfs/home/students/chit/is_benchmark/covid_config.sh


##### run down
# readfilesdir="/nfs/scratch/chit/isbench_covid/03_fastq_files_ordered/"
# for sample in $(ls ${readfilesdir} | grep "fastq.gz" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | sort | uniq); do
#     #echo SHUFFLING ${sample}...
#     #seqkit shuffle -j 40 ${readfilesdir}/${sample}_1.fastq.gz > ${readfilesdir}/shuffled_fastq/${sample}_1.fastq.shuffled.gz
#     if ! test -s /nfs/scratch/chit/isbench_covid/03_fastq_files_ordered/shuffled_fastq/${sample}_2.fastq.gz.shuffled; 
#         then
#         echo shuffling ${sample}
#         h
#         # pigz ${readfilesdir}/${sample}_1.fastq
#         # pigz ${readfilesdir}/${sample}_2.fastq
#         echo done ${sample}
#     fi
# done
