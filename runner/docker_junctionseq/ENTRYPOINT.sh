#!/bin/bash 

source /MOUNT/config.sh

R --version
echo "#########################"
echo "Running JunctionSeq from R"
echo "#########################"

Rscript /dependencies.R 

for x in $(cat /MOUNT/compare_list.txt); do
    con1=$(echo ${x} | cut -d ';' -f 1)
    con2=$(echo ${x} | cut -d ';' -f 2)

    Rscript /junctionseq.R /MOUNT $con1 $con2
    echo done $con1 $con2
done 