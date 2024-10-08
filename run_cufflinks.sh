#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/nfs/scratch/chit/cufflink_0.out
#SBATCH --error=/nfs/scratch/chit/cufflink_0.err
#SBATCH --job-name=cufflinks0
#SBATCH --mem-per-cpu=40G
#SBATCH --partition=exbio-cpu
#SBATCH --time=4-00:00:00


PATH=$PATH:/nfs/scratch/chit/cufflinks-2.2.1-gffpatch.Linux_x86_64/
PATH=$PATH:/nfs/scratch/chit/samtools-1.13
PATH=$PATH:/usr/bin/tophat2

set -euo pipefail

for arg in "$@"; do
  shift
  case "$arg" in
    "--help") set -- "$@" "-h" ;;
    "--config") set -- "$@" "-c" ;;
    "--name") set -- "$@" "-n" ;;

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
while getopts h:c:n: flag;
  do
    case "${flag}" in
      h) usage; exit 0;;
      c) config=${OPTARG};;
      n) name=${OPTARG};;
  esac
done

dir=$name
source ${config}

echo $name
echo running cufflinks
 
! test -d ${outputdir}/results/star_cufflinks && mkdir -p ${outputdir}/results/star_cufflinks
ls ${outputdir}/alignments | parallel --will-cite -j $nCores "  
    echo {}
    ! test -d ${outputdir}/results/star_cufflinks/{} && mkdir -p ${outputdir}/results/star_cufflinks/{}
    if ! test -s ${outputdir}/resuts/star_cufflinks/{}/isoforms.fpkm_tracking;
      then
      /nfs/scratch/chit/cufflinks-2.2.1.Linux_x86_64/cufflinks -g ${gtf} -o ${outputdir}/results/star_cufflinks/{} ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam -p 50
    fi
"


# ls ${outputdir}/alignments | parallel --will-cite -j $nCores "  
#   perl ${path}/runner/cufftrim.pl /nfs/scratch/chit/ref/ensembl_98/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai ${outputdir}/results/star_cufflinks/{}/transcripts.gtf > ${outputdir}/results/star_cufflinks/{}/transcripts_trim.gtf
# "

# ls ${outputdir}/alignments | parallel --will-cite -j $nCores "  
#   python ${path}/runner/cuffdedup.py --gtf ${outputdir}/results/star_cufflinks/{}/transcripts.gtf --output ${outputdir}/results/star_cufflinks/{}/transcripts_dedup.gtf
# "


echo running cuffmerge...
find ${outputdir}/results -type f | grep "transcripts.gtf" > ${outputdir}/results/assembly_list.txt

/nfs/scratch/chit/cufflinks-2.2.1.Linux_x86_64/cuffmerge -g ${gtf} -s ${fasta} -o ${outputdir}/results -p 12 ${outputdir}/results/assembly_list.txt


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

if ! test -s ${outputdir}/results/cuffdiff_results/isoforms.count_tracking; then 
  /nfs/scratch/chit/cufflinks-2.2.1.Linux_x86_64/cuffdiff -L g1,g2 -u ${outputdir}/results/merged.gtf -p 12 -b ${fasta} ${group1l:1} ${group2l:1} -o ${outputdir}/results/cuffdiff_results
fi 
# echo running cufflinks - bowtie
# ls ${outputdir}/alignments | grep "sample" | parallel --will-cite -j $nCores "  
#     echo {}
#     mkdir ${outputdir}/results/bowtie_cufflinks/{}
#     cufflinks -g ${gtf} -o ${outputdir}/results/bowtie_cufflinks/{} ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam -p 50
# "

# find results/cufflinks/ -type f | grep "transcripts.gtf" > assembly_list.txt
popd
echo DONE
