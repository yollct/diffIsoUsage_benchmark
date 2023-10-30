#!/bin/bash 

source /MOUNT/config.sh
#ls -R /opt/conda/envs/r3.3/

echo "#########################"
echo "Running JunctionSeq from R"
echo "#########################"

Rscript /dependencies.R 

for x in $(cat /MOUNT/compare_list.txt); do
    con1=$(echo ${x} | cut -d ';' -f 1)
    con2=$(echo ${x} | cut -d ';' -f 2)
    cores=$(echo ${x} | cut -d ';' -f 3)

    Rscript /junctionseq.R /MOUNT $con1 $con2 $cores
    echo done $con1 $con2
done 