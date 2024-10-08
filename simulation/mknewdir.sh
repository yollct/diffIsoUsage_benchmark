#!/bin/bash

 

set -euo pipefail
# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
    "--help") set -- "$@" "-h" ;;
    "--dir") set -- "$@" "-d" ;;
    "--type") set -- "$@" "-t" ;;
    "--rep") set -- "$@" "-r" ;;
    "--depth") set -- "$@" "-v" ;;
    "--noise") set -- "$@" "-n" ;;
 
    *)        set -- "$@" "$arg"
  esac
done


while getopts h:d:t:r:v:n: flag;
  do
    case "${flag}" in
      h) usage; exit 0;;
      d) dir=${OPTARG};;
      t) type=${OPTARG};;
      r) rep=${OPTARG};;
      v) depth=${OPTARG};;
      n) noise=${OPTARG};;
  esac
done

outputdir="/nfs/scratch/chit/new_simulations/"

rm -rf ${outputdir}/${dir}
mkdir ${outputdir}/${dir}

if [ $type == "single" ]; then
    echo "simulating single"
    cp ${outputdir}/single_temp/compare_list.txt ${outputdir}/${dir}  
    cp ${outputdir}/single_temp/meta.txt ${outputdir}/${dir}
    cp ${outputdir}/single_temp/config.sh ${outputdir}/${dir}
    cp ${outputdir}/pair_temp/JunctionSeq.flat.gff.gz ${outputdir}/${dir}

    sed 's/{{dirname}}/'${dir}'/g' ${outputdir}/single_temp/config.sh > ${outputdir}/${dir}/config.sh
    mkdir ${outputdir}/${dir}/fastq_sim
    mkdir ${outputdir}/${dir}/results
    mkdir ${outputdir}/${dir}/rsem_sim

    
    Rscript /nfs/scratch/chit/covid_ts/mod_isoform_result.R ${dir} ${rep} ${noise}
    
    
    bash /nfs/scratch/chit/covid_ts/rsem_sim.sh --depth ${depth} --type ${rep} --dir ${dir}

    Rscript /nfs/scratch/chit/new_simulations/fix_groundtruth.R ${dir}
fi

if [ $type == "pair" ]; then
    echo "simulating paired"
    cp ${outputdir}/single_temp/compare_list.txt ${outputdir}/${dir}  
    cp ${outputdir}/pair_temp/meta.txt ${outputdir}/${dir}
    cp ${outputdir}/pair_temp/config.sh ${outputdir}/${dir}
    cp ${outputdir}/pair_temp/JunctionSeq.flat.gff.gz ${outputdir}/${dir}

    sed 's/{{dirname}}/${dir}/g' ${outputdir}/pair_temp/config.sh > ${outputdir}/${dir}/config.sh
    mkdir ${outputdir}/${dir}/fastq_sim
    mkdir ${outputdir}/${dir}/results
    mkdir ${outputdir}/${dir}/rsem_sim

    
    Rscript /nfs/scratch/chit/isbench_covid/mod_isoform_result.R ${dir} ${rep} ${noise}
    
    bash /nfs/scratch/chit/isbench_covid/rsem_sim.sh --depth ${depth} --type ${rep} --dir ${dir} 
    Rscript /nfs/scratch/chit/new_simulations/fix_groundtruth.R ${dir}
fi

