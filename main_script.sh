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
shopt -s extglob

! test -d ${outputdir}/results && mkdir -p ${outputdir}/results || true

echo "Checking parameters..."
pushd ${outputdir}

### make sure meta group has two level
groups=$(awk 'BEGIN{FS=OFS=" "}{ print $2 }' ${meta} | grep -v "group" | sort | uniq)
groupl=""
for u in ${groups}; do 
  groupl="${groupl} ${u}"
done
groupl=( $groupl )
echo comparing condition ${groupl[@]}

# ###################### simulation if sim==TRUE ###############################
# if [ ! $sim == 'false' ]
# then

#   mkdir ${outputdir} || true
  
#   echo "Simulated reads..."
#   Rscript ${path}/scripts/sim_polyester.R $ngenes $depth $switches $rep $outdir

#   checkfasta=$(ls | grep "fasta" | wc -l)
#   echo "Done simulating reads. ${checkfasta} fasta files are in simulated_reads folder."
# fi

# echo $(find ${readfilesdir} -name "*fastq.gz" )

### prepare exon gtf

if ! test -f ${index}/this_gtf.exon;
  then
  echo "preparing exon annotation"
  python ${path}/runner/seqGSEA/prepare_exon_annotation_ensembl.py ${gtf} ${index}/this_gtf.exon
fi

if ! test -s ${index}/JunctionSeq.flat.gff.gz;
  then
  java -jar /nfs/home/students/chit/hartleys-QoRTs-099881f/QoRTs.jar makeFlatGff \
    --stranded \
    ${gtf} \
    ${index}/JunctionSeq.flat.gff.gz
fi
###check if there are bam files
# # ./ -name '*fastq.gz' -nowarn | grep "fastq.gz" | awk -F- '{ print $1 }' | sed 's/.\///g' | awk -F_ '{print $1}' | sort | uniq
# for i in $(ls ${readfilesdir} | grep "*fastq.gz\|*fasta.gz\|*fq.gz" );
# do echo "uncompressing $i"
# pigz -d ${readfilesdir}/${i}
# done
# wait



##################
#INDEXES
###################
###################### STAR index ###############################
if ! test -f ${index}/genomeParameters.txt
	then mkdir -p ${index} || true # to allow mkdir to fail gracefully, we add "|| true"
		echo "Indexing reference genome with STAR"
		/nfs/home/students/chit/STAR-2.7.9a/bin/Linux_x86_64/STAR \
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

############## RSEM index ###################
if ! test -f ${index}/asim.grp
  then mkdir -p ${index} || true
    echo "Building RSEM index"
    rsem-prepare-reference --gtf ${gtf} \
      ${fasta} ${index}/asim -p ${nCores}
  else echo "Reusing RSEM index"
fi 
##############################################

############## Salmon index ###################
# if ! test -f ${index}/salmon_index;
#   then mkdir -p ${index}/salmon_index || true
#     echo "Building Salmon index"
#     grep '^>' <(cat ${fasta}) | cut -d ' ' -f 1 > ${index}/decoys.txt
#     sed -i.bak -e 's/>//g' ${index}/decoys.txt
#     cat ${transcript_fasta} ${fasta} > ${index}/gentrome.fa.gz
#     salmon index -t ${index}/gentrome.fa.gz -d ${index}/decoys.txt -p 4 -i ${index}/salmon_index 
#   else echo "Reusing Salmon index"
# fi 
##############################################

################# kallisto index ###############
if ! test -f ${index}/kallisto-index
 then echo "Indexing reference genome with Kallisto"
 	kallisto index -i ${index}/kallisto-index $transcript_fasta
 else echo "Index found at index/kallisto-index, reusing it:"
fi

python ${path}/scripts/t2g.py -f ${transcript_fasta} -o ${index}/tx2g.gtf
###################################################

fastaq=$(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq" | tail -c 9)
if [ $(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq" | tail -c 9 | head -c 1) != "f" ]; then
  fastaq=$(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq" | tail -c 3)
fi
echo $fastaq

## Checking if fastq files are paired or single.
for i in $(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | sort | uniq);
	do fastqcount=$(ls ${readfilesdir} | grep "${i}" | grep -v "results" | wc -l)
	echo "Checking paires state of sample: ${readfilesdir}/${i}*  | Number of Fastqfiles: $fastqcount"
	echo "###################"
	
	pairstate="paired"
	if [ $fastqcount -eq 2 ] && [ $pairstate != "single" ] ; then
		if [ -z ${prevcount+set} ] ; then prevcount=${fastqcount}; fi
		if ! [ $prevcount -eq $fastqcount ] ; then echo "takes either unpaired or paired reads per session, but not mixed" && exit 1 ; fi
	  echo 'Fastq files found are paired' && pairstate="paired"
	else
  	echo 'Fastq files found are not paired' && pairstate="single"
	fi
done

# echo concatenating samples: {} if needed
# if ! test -f {}_1.fastq; 
#   then zcat {}_L001_R1_001.fastq.gz {}_L002_R1_001.fastq.gz {}_L003_R1_001.fastq.gz {}_L004_R1_001.fastq.gz > {}_1.fastq;
#   zcat {}_L001_R2_001.fastq.gz {}_L002_R2_001.fastq.gz {}_L003_R2_001.fastq.gz {}_L004_R2_001.fastq.gz > {}_2.fastq;
# fi

STAR --version

ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | cut -d "." -f 1 | sort | uniq | parallel --will-cite -j $nCores "
  echo running star...
  # paired data
  if [ $pairstate == "paired" ]; then
    echo Transcript quantification with Salmon for Sample: {}
    ! test -d ${outputdir}/salmon_out/{} && mkdir -p ${outputdir}/salmon_out/{} || true # to allow mkdir to fail gracefully, we add '|| true'
    if ! test -f ${outputdir}/salmon_out/{}/quant.sf
      then salmon quant -i ${index}/salmon_index -p 12 --validateMappings -l IU -o ${outputdir}/salmon_out/{} -1 ${readfilesdir}/{}_1.${fastaq} -2 ${readfilesdir}/{}_2.${fastaq} 2> ${outputdir}/salmon_out.stdout
    fi

    ! test -d ${outputdir}/kallisto_out/{} && mkdir -p ${outputdir}/kallisto_out/{}
    if ! test -s ${outputdir}/kallisto_out/{}/abundance.tsv; then
      kallisto quant -i ${index}/kallisto-index -o ${outputdir}/kallisto_out/{} --bias ${readfilesdir}/{}_1.${fastaq} ${readfilesdir}/{}_2.${fastaq} 2> ${outputdir}/kallisto.stdout

    fi

    echo Genome Alignment with STAR for Sample: {}
    ! test -d ${outputdir}/alignments/ && mkdir -p ${outputdir}/alignments/ || true # to allow mkdir to fail gracefully, we add '|| true'
    ! test -d ${outputdir}/alignments/{} && mkdir -p ${outputdir}/alignments/{} || true

    if ! test -s ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam
      then
      echo {} running STAR
      STAR --genomeDir ${index} \
        --runThreadN 12 \
        --readFilesIn  ${readfilesdir}/{}_1.${fastaq} ${readfilesdir}/{}_2.${fastaq} \
        --readFilesCommand cat \
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

    echo {} running RSEM
    if ! test -s ${outputdir}/rsem_out/{}/{}.isoforms.results; then
      rsem-calculate-expression --star --star-path /nfs/home/students/chit/STAR-2.7.8a/bin/Linux_x86_64 --paired-end ${readfilesdir}/{}_1.${fastaq} ${readfilesdir}/{}_2.${fastaq}  ${index}/asim ${outputdir}/rsem_out/{}/{}
      echo DONE for RSEM {}
    fi

    ! test -d ${outputdir}/qort && mkdir -p ${outputdir}/qort || true
    ! test -d ${outputdir}/qort/{} && mkdir -p ${outputdir}/qort/{} || true

    if ! test -s ${outputdir}/qort/{}/QC.spliceJunctionCounts.knownSplices.txt.gz; then
      echo {} running QoRT
      java -jar /nfs/home/students/chit/hartleys-QoRTs-099881f/QoRTs.jar QC \
                --stranded \
                ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam \
                ${gtf} \
                ${outputdir}/qort/{}/
      echo DONE for QoRT {}
    fi

    ! test -d ${outputdir}/results/htseq_exon && mkdir -p ${outputdir}/results/htseq_exon || true

    if ! test -s ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.sam; 
      then
      samtools view -h ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam > ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.sam
    fi 

    python ${path}/runner/seqGSEA/count_in_exons.py -p yes ${index}/this_gtf.exon ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.sam ${outputdir} {} ${meta}
    echo DONE HTSeq {}

  else 
    echo running for single-end read data
    echo Transcript quantification with Salmon for Sample: {}
    
    ! test -d ${outputdir}/salmon_out/{} && mkdir -p ${outputdir}/salmon_out/{} || true # to allow mkdir to fail gracefully, we add '|| true'
    if ! test -f ${outputdir}/salmon_out/{}/quant.sf
      then salmon quant -i ${index}/salmon_index -p 8 -l SF -o ${outputdir}/salmon_out/{} --validateMappings -r ${readfilesdir}/{}.${fastaq} 2> ${outputdir}/salmon_out.stdout
    fi

    ! test -d ${outputdir}/kallisto_out/{} && mkdir -p ${outputdir}/kallisto_out/{}
    if ! test -s ${outputdir}/kallisto_out/{}/abundance.tsv; then
      kallisto quant -i ${index}/kallisto-index -o ${outputdir}/kallisto_out/{} --bias --single -t 8 -l 66 -s 35 ${readfilesdir}/{}.${fastaq} 2> ${outputdir}/kallisto.stdout
    fi

    echo Genome Alignment with STAR for Sample: {}
    ! test -d ${outputdir}/alignments/ && mkdir -p ${outputdir}/alignments/ || true # to allow mkdir to fail gracefully, we add '|| true'
    ! test -d ${outputdir}/alignments/{} && mkdir -p ${outputdir}/alignments/{} || true

    if ! test -s ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam 
      then
      echo {} running STAR
      STAR --genomeDir ${index} \
        --runThreadN 12 \
        --readFilesIn  ${readfilesdir}/{}.${fastaq} \
        --readFilesCommand cat \
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
    else 
      echo STAR is done {}. skipped
    fi
    
    # ! test -d ${outputdir}/rsem_out && mkdir -p ${outputdir}/rsem_out || true
    # ! test -d ${outputdir}/rsem_out/{} && mkdir -p ${outputdir}/rsem_out/{} || true

    

    echo {} running QoRT
    ! test -d ${outputdir}/qort && mkdir -p ${outputdir}/qort || true
    ! test -d ${outputdir}/qort/{} && mkdir -p ${outputdir}/qort/{} || true

    if ! test -s ${outputdir}/qort/{}/QC.geneCounts.txt.gz; then
      echo {} running QoRT
      java -jar /nfs/home/students/chit/hartleys-QoRTs-099881f/QoRTs.jar QC \
                --stranded \
                --singleEnded \
                ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam \
                ${gtf} \
                ${outputdir}/qort/{}/
      echo DONE for QoRT {}
    else 
      echo QoRT is done {}. skipped
    fi


    ! test -d ${outputdir}/results/htseq_exon && mkdir -p ${outputdir}/results/htseq_exon || true

    if ! test -s ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.sam; 
      then
      samtools view -h ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam > ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.sam
    fi 

    python ${path}/runner/seqGSEA/count_in_exons.py ${index}/this_gtf.exon ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.sam ${outputdir} {} ${meta}
    echo DONE HTSeq {}
  
  fi
  "
 
  
  
  # if ! test -f ${outputdir}/rsem_out/{}/{}.genes.results
  #   then
  # rsem-calculate-expression --star --star-path /nfs/home/students/chit/STAR-2.7.8a/bin/Linux_x86_64 --star-gzipped-read-file --paired-end ${readfilesdir}/{}_1.${fastaq} ${readfilesdir}/{}_2.${fastaq}  ${index}/asim ${outputdir}/rsem_out/{}/{}

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



if ! test -f ${outputdir}/results/salmon_count.csv;
  then
  echo "Extracting salmon reads..."
  python ${path}/scripts/extract_salmon.py --dir ${outputdir}/salmon_out --outputfile ${outputdir}/results/ --pattern $pattern --type TPM --gtf ${index}/tx2g.gtf --meta $meta
  python ${path}/scripts/extract_kallisto.py --dir ${outputdir}/kallisto_out --outputfile ${outputdir}/results/ --pattern $pattern --type TPM --gtf ${index}/tx2g.gtf --meta $meta
  python ${path}/scripts/extract_rsem.py --dir ${outputdir}/rsem_out --outputfile ${outputdir}/results/ --pattern $pattern --type TPM --gtf $gtf  --meta $meta
  echo "Finished running salmon."
fi



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
if ! test -s ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! test -s ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt || ! test -s ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt; then
  Rscript ${path}/runner/drimseq.R $outputdir $path $meta ${groupl[@]}
fi

echo "running iso-KTSP"
if ! test -f ${outputdir}/results/salmon_isoktsp_output_${groupl[0]}_${groupl[1]}.txt; then
  java -jar /nfs/home/students/chit/iso-kTSP_v1.0.3.jar ${outputdir}/results/salmon_isoktsp_count.txt -i -n 2 -s 5000 -o ${outputdir}/results/salmon_isoktsp_output_${groupl[0]}_${groupl[1]}.txt --seed 1234 -c ${groupl[0]} ${groupl[1]}
fi

if ! test -f ${outputdir}/results/kal_isoktsp_output_${groupl[0]}_${groupl[1]}.txt;
  then
  java -jar /nfs/home/students/chit/iso-kTSP_v1.0.3.jar ${outputdir}/results/kal_isoktsp_count.txt -i -n 2 -s 5000 -o ${outputdir}/results/kal_isoktsp_output_${groupl[0]}_${groupl[1]}.txt --seed 1234 -c ${groupl[0]} ${groupl[1]}
fi

echo "Finished running iso-KTSP"
echo "#########################"
echo "Running DEXSeq from R"
echo "#########################"
if ! grep -q "dexseq" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "dexseq" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "dexseq" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/dexseq.R $outputdir $path $meta ${groupl[@]}
fi

echo "#########################"
echo "Running DTurtle from R"
echo "#########################"
if ! grep -q "dturtle" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "dturtle" ${outputdir}/results/kal_res_tx_${groupl[0]}_${groupl[1]}.txt || ! grep -q "dturtle" ${outputdir}/results/rsem_res_tx_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/dtu.R $outputdir $path $meta $gtf ${groupl[@]}
fi
echo "#########################"
echo "Running iso-KTSP from R"
echo "#########################"
if ! grep -q "iso_ktsp" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "iso_ktsp" ${outputdir}/results/kal_res_tx_${groupl[0]}_${groupl[1]}.txt || ! grep -q "iso_ktsp" ${outputdir}/results/rsem_res_tx_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/iso-ktsp.R $outputdir $path ${groupl[@]}
fi


#Rscript cuffdiff.R

echo "#########################"
echo "Running seqGSEA from R"
echo "#########################"
if ! grep -q "seqGSEA" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "seqGSEA" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "seqGSEA" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; then
  for file in $(ls ${outputdir}/results/htseq_exon); 
  do
    sed 's/"//g' ${outputdir}/results/htseq_exon/${file} > ${outputdir}/results/htseq_exon/temp.txt 
    mv ${outputdir}/results/htseq_exon/temp.txt ${outputdir}/results/htseq_exon/${file}
  done
  Rscript ${path}/runner/seqGSEA/seqgsea.R $outputdir $path $meta $index ${groupl[@]}
fi

echo "#########################"
echo "Running JunctionSeq from R"
echo "#########################"
if ! grep -q "junctionseq" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "junctionseq" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "junctionseq" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
  then
  #docker run -v ${outputdir}:/MOUNT --rm --name 'jcseq' jcseq 
  Rscript ${path}/runner/junctionseq_res.R $outputdir $path $meta ${groupl[@]} 
fi

echo "#########################"
echo "Running Cuffdiff from R"
echo "#########################"
if ! grep -q "cuffdiff" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "cuffdiff" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "cuffdiff" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/cuffdiff.R $outputdir $path $meta ${groupl[@]} $gtf
fi


echo "#########################"
echo "Running Saturn from R"
echo "#########################"
if ! grep -q "saturn" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "saturn" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "saturn" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/saturn.R $outputdir $path $meta ${groupl[@]} 
fi