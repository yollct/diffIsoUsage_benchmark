#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output=isodot.err
#SBATCH --error=isodot.out
#SBATCH --job-name=isodot
#SBATCH --array=1-10
#SBATCH --mem-per-cpu=20G
#SBATCH --time=48:00:00

path="/nfs/home/students/chit/is_benchmark/results"
output="/nfs/home/students/chit/is_benchmark/"
PATH=$PATH:/nfs/home/students/chit/samtools-1.13

cd ${output}/alignments
# sample="sample_01"
for sample in $(ls | grep "sample");
do 
    #cd $sample

# #samtools sort -n -o ${sample}Aligned.sortedByCoord.out.sorted.bam -@ 10 ${sample}Aligned.sortedByCoord.out.bam
# if [ ! -f ${sample}Aligned.sortedByCoord.out.bam.bai ]; then
#     echo "indexing bam file..."
#     samtools index -@ 20 ${sample}Aligned.sortedByCoord.out.bam 
# fi

# #samtools sort -n -o ${sample}Aligned.sortedByCoord.out.sorted.bam -@ 10 ${sample}Aligned.sortedByCoord.out.bam
# bam="${output}/alignments/${sample}/${sample}Aligned.sortedByCoord.out.sorted.bam"
    ##make frag size
    pushd ${output}
    if [ ! -f ./alignments/${sample}/${sample}_isodot_output_knowniso.RData ];then
        # echo "making frag size file"
        # samtools view -@ 10 ./alignments/${sample}/${sample}Aligned.sortedByCoord.out.bam  \
        # | awk '{ print ($8 >= $4) ? $8-$4+76 : $4-$8+76 }' \
        # | sort -n | uniq -c > ./alignments/${sample}/${sample}_insertSize_dist.txt
        # ls ${sample} -lah

        Rscript ${output}/runner/isodot/isodot.R ${sample}
        echo "done for ${sample}"
    else 
        echo "${sample} is already done"
    fi
    
    popd

    
#samtools view  -@ 10 ${sample}Aligned.sortedByCoord.out.sorted.bam > ${output}/alignments/${sample}/${sample}_insertSize_dist.txt
# cd ..
done

Rscript ${output}/runner/isodot/isodot_dtu.R 

