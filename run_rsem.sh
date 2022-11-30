#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=/nfs/scratch/chit/rsem.out
#SBATCH --error=/nfs/scratch/chit/rsem.err
#SBATCH --job-name=rsem
#SBATCH --mem-per-cpu=20G
#SBATCH --time=14-00:00:00


PATH=$PATH:/nfs/home/students/chit/RSEM-1.3.2/

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
 
! test -d ${outputdir}/rsem_out && mkdir -p ${outputdir}/rsem_out

fastaq=$(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq" | tail -c 9)
if [ $(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq" | tail -c 9 | head -c 1) != "f" ]; then
  fastaq=$(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq" | tail -c 3)
fi
echo $fastaq

ls ${outputdir}/alignments | grep "SRR" | parallel --will-cite -j $nCores "  
    echo {}
    ! test -d ${outputdir}/rsem_out/{} && mkdir -p ${outputdir}/rsem_out/{}
    echo {} running RSEM
    if ! test -s ${outputdir}/rsem_out/{}/{}.isoforms.results; then
      rsem-calculate-expression --star --star-path /nfs/home/students/chit/STAR-2.7.8a/bin/Linux_x86_64  ${readfilesdir}/{}.${fastaq} ${index}/asim ${outputdir}/rsem_out/{}/{}
      echo DONE for RSEM {}
    else 
      echo RSEM is done {}. skipped
    fi
"
