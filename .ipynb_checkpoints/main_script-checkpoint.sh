#!/bin/bash
PATH=$PATH:/nfs/home/students/chit/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/nfs/home/students/chit/samtools-1.13
PATH=$PATH:/nfs/home/students/chit/RSEM-1.3.2/
PATH=$PATH:/nfs/home/students/chit/tophat-2.1.1.Linux_x86_64
PATH=$PATH:/nfs/home/students/chit/bowtie2-2.4.5-linux-x86_64
PATH=$PATH:/nfs/home/students/chit/STAR-2.7.8a/bin/Linux_x86_64

set -euo pipefail
# Transform long options to short ones
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

echo "$config"
source ${config}

! test -d ${outputdir}/results && mkdir -p ${outputdir}/results || true
batches=$((rep*2))
echo $batches

echo "Checking parameters..."
pushd ${outputdir}

###################### simulation if sim==TRUE ###############################
if [ ! $sim == 'false' ]
then

  mkdir ${outputdir} || true
  
  echo "Simulated reads..."
  Rscript ${path}/scripts/sim_polyester.R $ngenes $depth $switches $rep $outdir

  checkfasta=$(ls | grep "fasta" | wc -l)
  echo "Done simulating reads. ${checkfasta} fasta files are in simulated_reads folder."
fi

echo $(find ${readfilesdir} -name "*fastq.gz" )

###check if there are bam files
#find ./ -name '*fastq.gz' -nowarn | grep "fastq.gz" | awk -F- '{ print $1 }' | sed 's/.\///g' | awk -F_ '{print $1}' | sort | uniq
# for i in $(find ${readfilesdir} -name "*fastq.gz" );
# do echo "uncompressing $i"
# pigz -d $i
# done
# wait

##################
#INDEXES
###################
###################### STAR index ###############################
if ! test -f ${index}/genomeParameters.txt
	then mkdir -p ${index} || true # to allow mkdir to fail gracefully, we add "|| true"
		echo "Indexing reference genome with STAR"
		STAR \
		--runMode genomeGenerate \
		--genomeDir ${index} \
		--genomeFastaFiles $fasta \
		--runThreadN $nCores \
		--sjdbGTFfile $gtf \
		--outFileNamePrefix ${outputdir}/alignments
	else echo "Indexes found at index/starindex, reusing them:"
    ls ${index} -lah
fi 

################# Build bowtie index #######################
# if ! test -d ${bowtie_index}/genome_index
# 	then mkdir -p ${bowtie_index} || true # to allow mkdir to fail gracefully, we add "|| true"
#     mkdir -p ${bowtie_index}/genome_index || true
# 		echo "Indexing reference genome with Bowtie"
# 		bowtie2-build \
#       --threads $nCores \
#       ${fasta} \
#       ${bowtie_index}/genome_index/genome_in
# 	else echo "Indexes found at index/bowtie_index, reusing them:"
# fi

# if ! test -d ${bowtie_index}/transcript_index
# 	then mkdir -p ${bowtie_index} || true # to allow mkdir to fail gracefully, we add "|| true"
#     mkdir -p ${bowtie_index}/transcript_index || true
# 		echo "Indexing reference transcriptome with Bowtie"
# 		bowtie2-build \
#       --threads $nCores \
#       ${transcript_fasta} \
#       ${bowtie_index}/transcript_index/trans_in
# 	else echo "Indexes found at index/bowtie_index, reusing them:"
# fi
#########################################################
echo hi
############## RSEM index ###################
if ! test -f ${index}/asim.grp
  then mkdir -p ${index} || true
    echo "Building RSEM index"
    rsem-prepare-reference --gtf ${gtf} \
      ${fasta} ${index}/asim -p ${nCores}
  else echo "Reusing RSEM index"
fi 
##############################################

# echo concatenating samples: {} if needed
# if ! test -f {}_1.fastq; 
#   then zcat {}_L001_R1_001.fastq.gz {}_L002_R1_001.fastq.gz {}_L003_R1_001.fastq.gz {}_L004_R1_001.fastq.gz > {}_1.fastq;
#   zcat {}_L001_R2_001.fastq.gz {}_L002_R2_001.fastq.gz {}_L003_R2_001.fastq.gz {}_L004_R2_001.fastq.gz > {}_2.fastq;
# fi

fastaq=$(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz" | tail -c 9)
echo $fastaq
ls ${readfilesdir} | grep "fasta.gz\|fastq.gz" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | sort | uniq | parallel --will-cite -j $nCores "
  echo running star...

	echo Transcript quantification with Salmon for Sample: {}
	! test -d ${outputdir}/salmon_out/{} && mkdir -p ${outputdir}/salmon_out/{} || true # to allow mkdir to fail gracefully, we add '|| true'
	if ! test -f ${outputdir}/salmon_out/{}/quant.sf
		then salmon quant -i ${outputdir}/salmon_index -o ${outputdir}/salmon_out/{} --bias ${readfilesdir}/{}_1.fastq.gz ${readfilesdir}/{}_2.fastq.gz 2> ${outputdir}/salmon_out.stdout
	fi

	echo Genome Alignment with STAR for Sample: {}
	! test -d ${outputdir}/alignments/ && mkdir -p ${outputdir}/alignments/ || true # to allow mkdir to fail gracefully, we add '|| true'
  ! test -d ${outputdir}/alignments/{} && mkdir -p ${outputdir}/alignments/{} || true



  if ! test -s ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam
    then
    echo {} running STAR
    STAR --genomeDir ${index} \
      --runThreadN 10 \
      --readFilesIn  ${readfilesdir}/{}_1.${fastaq} ${readfilesdir}/{}_2.${fastaq} \
      --readFilesCommand zcat \
      --outSAMstrandField intronMotif \
      --outFilterIntronMotifs RemoveNoncanonical \
      --outSAMattrIHstart 0 \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes jM jI\
      --alignSoftClipAtReferenceEnds No \
      --quantMode TranscriptomeSAM \
      --sjdbGTFfile $gtf \
      --outFileNamePrefix ${outputdir}/alignments/{}/{}

    echo DONE STAR for {}
  fi
	
	
  ! test -d ${outputdir}/rsem_out && mkdir -p ${outputdir}/rsem_out || true
  ! test -d ${outputdir}/rsem_out/{} && mkdir -p ${outputdir}/rsem_out/{} || true
  if ! test -f ${outputdir}/rsem_out/{}/{}.genes.results
    then
    echo {} running RSEM
    rsem-calculate-expression --star --star-path /nfs/home/students/chit/STAR-2.7.8a/bin/Linux_x86_64 --star-gzipped-read-file --paired-end ${readfilesdir}/{}_1.${fastaq} ${readfilesdir}/{}_2.${fastaq}  ${index}/asim ${outputdir}/rsem_out/{}/{}
    echo DONE for RSEM {}
  fi
	"

  #${outputdir}/alignments/{}/{}Aligned.toTranscriptome.out.bam
  ### tophat code
  #   echo Alignment with Tophat for Sample: {}
  # ! test -d ${outputdir}/alignments/tophat && mkdir -p ${outputdir}/alignments/tophat || true
  # ! test -d ${outputdir}/alignments/tophat/{} && mkdir -p ${outputdir}/alignments/tophat/{} || true


  # if ! test -f ${readfilesdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam.bai
	#   then echo indexing {}
	# 	samtools index ${readfilesdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam ${readfilesdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam.bai
	# fi
echo hi

if ! test -d ${outputdir}/results/salmon_count.csv
  then
  echo "Extracting salmon reads..."
  python ${path}/scripts/extract_salmon.py --dir ${outputdir}/salmon_out --outputfile ${outputdir}/results/ --pattern $pattern --type TPM --gtf $gtf --meta $meta
  echo "Finished running salmon."
fi

echo "running iso-KTSP"
# if ! test -d ${outputdir}/results/isoktsp_output.txt
#   then
#   #java -jar /nfs/home/students/chit/iso-kTSP_v1.0.3.jar ${outputdir}/results/isoktsp_count.txt -i -n 2 -s 5000 -o ${outputdir}/results/isoktsp_output.txt --seed 1234 -c "1" "2"
# fi
echo "Finished running iso-KTSP"

#python ${path}/scripts/extract_salmon.py --dir ${readfilesdir}/salmon_out --outputfile ${readfilesdir}/results/salmon_count.csv --pattern sample --type TPM


echo "run isodot"

# ls ${readfilesdir}/salmon_out |  parallel --will-cite -j $nCores "

#   if ! test -f ${readfilesdir}/alignments/{}/{}_isodot_output_knowniso.RData 
#     then echo making frag size file: {}
#     samtools view -@ 10 ${readfilesdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam  \
#     | awk '{ print ($8 >= $4) ? $8-$4+76 : $4-$8+76 }' \
#     | sort -n | uniq -c > ${readfilesdir}/alignments/{}/{}_insertSize_dist.txt

#     Rscript ${path}/runner/isodot/isodot.R {}
#     echo done for {}
#   else 
#     echo {} is already done
#   fi"

# if [ -f ${readfilesdir}/results/cufflinks/{}/isoform.fpkm_tracking ]; then
#         echo {} exists
#         continue
#     else

###################### Run cufflinks for STAR resultss #######
#run cufflink in separate scripts 'cos it takes too long.....
############################################################


echo "Done preprocessing."
echo $sim
if [ $sim == 'true' ];then
  Rscript ${path}/runner/groundtruth.R $outputdir
fi

echo "#########################"
echo "Running DRIMSeq from R"
echo "#########################"
Rscript ${path}/runner/drimseq.R $outputdir $path $meta

echo "#########################"
echo "Running DEXSeq from R"
echo "#########################"
Rscript ${path}/runner/dexseq.R $outputdir $path $meta 

echo "#########################"
echo "Running DTurtle from R"
echo "#########################"
Rscript ${path}/runner/dtu.R $outputdir $path $meta $gtf

echo "#########################"
echo "Running iso-KTSP from R"
echo "#########################"
Rscript ${path}/runner/iso-ktsp.R $outputdir
#Rscript cuffdiff.R


popd