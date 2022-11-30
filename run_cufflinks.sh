#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=/nfs/scratch/chit/cufflink_s2.out
#SBATCH --error=/nfs/scratch/chit/cufflink_s2.err
#SBATCH --job-name=s2_cufflink
#SBATCH --mem-per-cpu=20G
#SBATCH --time=14-00:00:00


PATH=$PATH:/nfs/home/students/chit/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/nfs/home/students/chit/samtools-1.13
PATH=$PATH:/nfs/home/students/chit/RSEM/bin/
PATH=$PATH:/usr/bin/tophat2

set -euo pipefail

for arg in "$@"; do
  shift
  case "$arg" in
    "--help") set -- "$@" "-h" ;;
    "--config") set -- "$@" "-c" ;;

    *)        set -- "$@" "$arg"
  esac
done

usage(){
  echo "ISbenchmark - run simulations and benchmark isoform switch tools"
  echo " "
  echo "[options] application [arguments]"
  echo " "
  echo "options:"
  echo "-h,         show brief help"
  echo "-c,         config file path"
}
exit_abnormal() {                         # Function: Exit with error.
  usage
  exit 1
}

      ######## simulation ########
while getopts h:c: flag;
  do
    case "${flag}" in
      h) usage; exit 0;;
      c) config=${OPTARG};;
  esac
done
source ${config}

echo running cufflinks
 
! test -d ${outputdir}/results/star_cufflinks && mkdir -p ${outputdir}/results/star_cufflinks

ls ${outputdir}/alignments | grep "SRR" | parallel --will-cite -j $nCores "  
    echo {}
    ! test -d ${outputdir}/results/star_cufflinks/{} && mkdir -p ${outputdir}/results/star_cufflinks/{}
    if ! test -s ${outputdir}/resuts/star_cufflinks/{}/isoforms.fpkm_tracking;
      then
      cufflinks -g ${gtf} -o ${outputdir}/results/star_cufflinks/{} ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam -p 50
    fi
"


echo running cuffmerge...
find ${outputdir}/results -type f | grep "transcripts.gtf" > ${outputdir}/results/assembly_list.txt

cuffmerge -g ${gtf} -s ${fasta} -o ${outputdir}/results -p 4 ${outputdir}/results/assembly_list.txt


echo here
pushd ${outputdir}
echo running cuffdiff...
######### change here to fit the group name of meta data ###########
group1=$(awk 'BEGIN{FS=OFS=" "}{ if ($2 == "T") print"/alignments/"$1"/"$1"Aligned.sortedByCoord.out.bam" }' ${meta} | grep -v "group") 
group2=$(awk 'BEGIN{FS=OFS=" "}{ if ($2 == "N") print"/alignments/"$1"/"$1"Aligned.sortedByCoord.out.bam" }' ${meta} | grep -v "group") 
group1l=""

for i in $group1;
do 
    group1l="${group1l},${outputdir}${i}"
done;

group2l=""
for i in $group2;
do 
    group2l="${group2l},${outputdir}${i}"
done;

echo ${group1l:1}
cuffdiff -L g1,g2 -u ${outputdir}/results/merged.gtf -p 4 -b ${fasta} ${group1l:1} ${group2l:1} -o ${outputdir}/results/cuffdiff_results
# echo running cufflinks - bowtie
# ls ${outputdir}/alignments | grep "sample" | parallel --will-cite -j $nCores "  
#     echo {}
#     mkdir ${outputdir}/results/bowtie_cufflinks/{}
#     cufflinks -g ${gtf} -o ${outputdir}/results/bowtie_cufflinks/{} ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam -p 50
# "

# find results/cufflinks/ -type f | grep "transcripts.gtf" > assembly_list.txt
popd
echo DONE