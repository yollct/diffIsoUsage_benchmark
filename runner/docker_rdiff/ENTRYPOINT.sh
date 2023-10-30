#!/bin/bash 

source /MOUNT/config.sh

echo "#########################"
echo "Running rDiff"
echo "#########################"

# Rscript /dependencies.R 

# for x in $(cat /MOUNT/compare_list.txt); do
#     con1=$(echo ${x} | cut -d ';' -f 1)
#     con2=$(echo ${x} | cut -d ';' -f 2)

#     Rscript /iuta.R /MOUNT $con1 $con2
#     echo done $con1 $con2
# done 
ls -R /usr/lib/cruft 

bams='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
samples=$(ls /MOUNT/fastq_sim | grep "fasta.gz\|fastq.gz\|fq\|fastq" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | cut -d "." -f 1 | sort | uniq)

mkdir /MOUNT/results/rdiff/ 

group1=$(awk 'BEGIN{FS=OFS=" "}{ if ($2 == "T") print"/MOUNT/alignments/"$1"/"$1"Aligned.sortedByCoord.out.bam" }' /MOUNT/meta.txt | grep -v "group") 
group2=$(awk 'BEGIN{FS=OFS=" "}{ if ($2 == "N") print"/MOUNT/alignments/"$1"/"$1"Aligned.sortedByCoord.out.bam" }' /MOUNT/meta.txt | grep -v "group") 

for i in $group1;
do 
    group1l="${group1l},${i}"
done;

echo $group1l
group2l=""
for i in $group2;
do 
    group2l="${group2l},${i}"
done;

which octave 

/rDiff/bin/rdiff -o /MOUNT/results/ -a ${group1l:1} -b ${group2l:1} -g /ref/Homo_sapiens.GRCh38.110.gff3 -m param -L 100 -m 30

