#!/bin/bash 

source /MOUNT/config.sh


R --version
echo "#########################"
echo "Running rSeqDiff from R"
echo "#########################"

#Rscript /dependencies.R 

#/rseq-0.2.2-src/rseq gen_exons_junctions -v 2 /ref/Homo_sapiens.GRCh38.98.refflat /ref/chrlist.txt

mkdir /MOUNT/rseqdiff/

# ls /MOUNT/fastq_sim | grep "fasta.gz\|fastq.gz\|fq" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | cut -d "." -f 1 | sort | uniq | parallel --will-cite -j 8 "
#     echo running seq map {}
#     /rseq-0.2.2-src/seqmap 2 /MOUNT/fastq_sim/{}.fq /ref/Homo_sapiens.GRCh38.98.refflat.exons_junctions.fa /MOUNT/rseqdiff/{}_out.txt /eland:3
#     #/rseq-0.2.2-src/seqmap 2 /MOUNT/fastq_sim/{}.1_fq /ref/Homo_sapiens.GRCh38.98.refflat.exons_junctions.fa /MOUNT/rseqdiff/{}_out2.txt /eland:3
#     echo done seqmap {}

#     echo running rseq {}
#     /rseq-0.2.2-src/rseq comp_exp -v 2 /ref/Homo_sapiens.GRCh38.98.refflat.subexons.txt /MOUNT/rseqdiff/{}_out.txt 
#     echo done rseq {}
# "

for x in $(cat /MOUNT/compare_list.txt); do
    con1=$(echo ${x} | cut -d ';' -f 1)
    con2=$(echo ${x} | cut -d ';' -f 2)

    Rscript /rseqdiff.R /MOUNT $con1 $con2
    echo done $con1 $con2
done 