#!/bin/bash

PATH=$PATH:/opt/anaconda3/bin/
PATH=$PATH:/nfs/scratch/chit/cufflinks-2.2.1.Linux_x86_64
PATH=$PATH:/nfs/scratch/chit/samtools-1.13
PATH=$PATH:/nfs/scratch/chit/RSEM-1.3.2/
PATH=$PATH:/nfs/scratch/chit/tophat-2.1.1.Linux_x86_64
PATH=$PATH:/nfs/scratch/chit/bowtie2-2.4.5-linux-x86_64
PATH=$PATH:/nfs/scratch/chit/STAR-2.7.9a/bin/Linux_x86_64
PATH=$PATH:/nfs/scratch/chit/DSGseq/DSG-0.1.0/SeqExpress
PATH=$PATH:/nfs/scratch/chit/bedtools2/bin
#PATH=$PATH:/nfs/home/students/chit/.conda/envs/nease/pkgs/r-base-3.6.1-haffb61f_2/bin

set -euo pipefail
# Transform long options to short ones
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

echo "$config"
echo $name
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
  java -jar /nfs/scratch/chit/hartleys-QoRTs-099881f/QoRTs.jar makeFlatGff \
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
		/nfs/scratch/chit/STAR-2.7.9a/bin/Linux_x86_64/STAR \
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

# # ACTIVATE ANACONDA

# echo "conda activated"

# chrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22'
# for thistmp in $(find ${outputdir}/alignments/ | grep "Aligned.sortedByCoord.out.bam$" | grep -v "tmp"); do
#   for chr in ${chrs}; do
#     echo separate bam $chr

#     if ! test -s ${thistmp}.bai; then
#       samtools index -@ 12 $thistmp 
#     fi 
#     samtools view -b -h -o ${thistmp}_tmp${chr}.bam ${thistmp} $chr
#     # (samtools view -H ${thistmp}; samtools view ${thistmp} | awk -v a=$chr '/^@/ || $3 == a ') | samtools view -bo ${thistmp}_tmp${chr}.bam &
#   done 
# done

# for thistmp in $(find ${outputdir}/alignments/ | grep "Aligned.sortedByCoord.out.bam$"); do
#   (samtools view -H ${thistmp}; samtools view ${thistmp} | awk '/^@/ || !($3 >= 1 && $3 <23)') | samtools view -bo ${thistmp}_tmp23.bam
# done

bams='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23'
samples=$(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq\|fastq" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | cut -d "." -f 1 | sort | uniq)


parallel --will-cite -j $nCores "
  echo running star...

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
        --runThreadN 8 \
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
        --twopassMode Basic \
        --outFileNamePrefix ${outputdir}/alignments/{}/{}

      echo DONE STAR for {}
    fi

    # if ! test -s ${outputdir}/alignments/{}/{}Aligned.toTranscriptome.out.bam.bai; then
    #   samtools index -@ 8 ${outputdir}/alignments/{}/{}Aligned.toTranscriptome.out.bam
    # fi

    ! test -d ${outputdir}/rsem_out && mkdir -p ${outputdir}/rsem_out || true
    ! test -d ${outputdir}/rsem_out/{} && mkdir -p ${outputdir}/rsem_out/{} || true

    echo {} running RSEM
    if ! test -s ${outputdir}/rsem_out/{}/{}.isoforms.results; then
      /nfs/scratch/chit/RSEM-1.3.3/rsem-calculate-expression --star --star-path /nfs/scratch/chit/STAR-2.7.9a/bin/Linux_x86_64 --paired-end ${readfilesdir}/{}_1.${fastaq} ${readfilesdir}/{}_2.${fastaq}  ${index}/asim ${outputdir}/rsem_out/{}/{}
      echo DONE for RSEM {}
    fi

  else 
    echo running for single-end read data
    echo Transcript quantification with Salmon for Sample: {1}
    
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
    rm -rf ${outputdir}/alignments/{}/tmp/ && mkdir ${outputdir}/alignments/{}/tmp/
    chmod -R 777 ${outputdir}/alignments/{}/tmp/

    if ! test -s ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam 
      then
      echo {} running STAR
      STAR --genomeDir ${index} \
        --runThreadN 4 \
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
        --twopassMode Basic \
        --outFileNamePrefix ${outputdir}/alignments/{}/{} 

      echo DONE STAR for {}
    else 
      echo STAR is done {}. skipped
    fi

   

    ! test -d ${outputdir}/rsem_out && mkdir -p ${outputdir}/rsem_out || true
    ! test -d ${outputdir}/rsem_out/{} && mkdir -p ${outputdir}/rsem_out/{} || true
    
    if ! test -f ${outputdir}/rsem_out/{}/{}.genes.results
    then
      /nfs/scratch/chit/RSEM-1.3.3/rsem-calculate-expression --star --star-path /nfs/scratch/chit/STAR-2.7.9a/bin/Linux_x86_64 ${readfilesdir}/{}.${fastaq} ${index}/asim ${outputdir}/rsem_out/{}/{}
  
    fi
      
  fi
  " ::: $samples 


parallel --will-cite -j $nCores "
  echo running qort
  if [ $pairstate == "paired" ]; then
    ! test -d ${outputdir}/qort && mkdir -p ${outputdir}/qort || true
    ! test -d ${outputdir}/qort/{1} && mkdir -p ${outputdir}/qort/{1} || true

    if ! test -s ${outputdir}/qort/{}/QC.spliceJunctionCounts.knownSplices.txt.gz; then
      echo {} running QoRT
      java -jar /nfs/scratch/chit/hartleys-QoRTs-099881f/QoRTs.jar QC \
                --stranded \
                ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam \
                ${gtf} \
                ${outputdir}/qort/{}/
      echo DONE for QoRT {}
    fi
    
    # echo {1} running QoRT sub {2}
    # ! test -d ${outputdir}/qort{2} && mkdir -p ${outputdir}/qort{2} || true
    # ! test -d ${outputdir}/qort{2}/{1} && mkdir -p ${outputdir}/qort{2}/{1} || true


    # if ! test -s ${outputdir}/qort{2}/{1}/QC.geneCounts.txt.gz; then
    #   echo ${outputdir}/alignments/{1}/{1}Aligned.sortedByCoord.out.bam_tmp{2}.bam
    #   java -jar /nfs/scratch/chit/hartleys-QoRTs-099881f/QoRTs.jar QC \
    #             --stranded \
    #             ${outputdir}/alignments/{1}/{1}Aligned.sortedByCoord.out.bam_tmp{2}.bam \
    #             ${gtf} \
    #             ${outputdir}/qort{2}/{1}/
    #   echo DONE for QoRT {1}
    # else 
    #   echo QoRT is done {1}. skipped
    # fi
  else 
     # if ! test -s ${outputdir}/alignments/{1}/{1}Aligned.toTranscriptome.out.bam.bai; then
    #   samtools sort -@ 8 -o ${outputdir}/alignments/{1}/{1}Aligned.toTranscriptome.sorted.out.bam ${outputdir}/alignments/{1}/{1}Aligned.toTranscriptome.out.bam
    #   samtools index -@ 8 ${outputdir}/alignments/{1}/{1}Aligned.toTranscriptome.sorted.out.bam
    # fi
    
    echo {} running QoRT
    ! test -d ${outputdir}/qort && mkdir -p ${outputdir}/qort || true
    ! test -d ${outputdir}/qort/{1} && mkdir -p ${outputdir}/qort/{} || true

    if ! test -s ${outputdir}/qort/{}/QC.geneCounts.txt.gz; then
      echo {} running QoRT
      java -jar /nfs/scratch/chit/hartleys-QoRTs-099881f/QoRTs.jar QC \
                --stranded \
                --singleEnded \
                ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam \
                ${gtf} \
                ${outputdir}/qort/{}/
      echo DONE for QoRT {}
    else 
      echo QoRT is done {}. skipped
    fi

    # echo {} running QoRT sub 1
    # ! test -d ${outputdir}/qort1 && mkdir -p ${outputdir}/qort1 || true
    # ! test -d ${outputdir}/qort1/{} && mkdir -p ${outputdir}/qort1/{} || true


    # if ! test -s ${outputdir}/qort1/{}/QC.geneCounts.txt.gz; then
    #   echo {} running QoRT
    #   java -jar /nfs/scratch/chit/hartleys-QoRTs-099881f/QoRTs.jar QC \
    #             --stranded \
    #             --singleEnded \
    #             ${outputdir}/alignments/{}/{}Aligned.sortedByCoord.out.bam_tmp1.bam \
    #             ${gtf} \
    #             ${outputdir}/qort1/{}/
    #   echo DONE for QoRT {}
    # else 
    #   echo QoRT is done {}. skipped
    # fi

    # echo {1} running QoRT sub {2}
    # ! test -d ${outputdir}/qort{2} && mkdir -p ${outputdir}/qort{2} || true
    # ! test -d ${outputdir}/qort{2}/{1} && mkdir -p ${outputdir}/qort{2}/{1} || true


    # echo {1} running QoRT sub {2}
    # ! test -d ${outputdir}/qort{2} && mkdir -p ${outputdir}/qort{2} || true
    # ! test -d ${outputdir}/qort{2}/{1} && mkdir -p ${outputdir}/qort{2}/{1} || true


    # if ! test -s ${outputdir}/qort{2}/{1}/QC.geneCounts.txt.gz; then
    #   echo ${outputdir}/alignments/{1}/{1}Aligned.sortedByCoord.out.bam_tmp{2}.bam
    #   java -jar /nfs/scratch/chit/hartleys-QoRTs-099881f/QoRTs.jar QC \
    #             --stranded \
    #             --singleEnded \
    #             ${outputdir}/alignments/{1}/{1}Aligned.sortedByCoord.out.bam_tmp{2}.bam \
    #             ${gtf} \
    #             ${outputdir}/qort{2}/{1}/
    #   echo DONE for QoRT {1}
    # else 
    #   echo QoRT is done {1}. skipped
    # fi

  fi

" ::: $samples 


for fas in $(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | cut -d "." -f 1 | sort | uniq); do
  echo running DSGseq...

  ! test -d ${outputdir}/results/dsgseq && mkdir -p ${outputdir}/results/dsgseq || true

  if ! test -s ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.bed; 
    then
    bamToBed -i ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.bam > ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.bed &
  fi 
  
  if ! test -s ${outputdir}/results/dsgseq/${fas}.count; then
    SeqExpress count ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.bed ${genepred} ${outputdir}/results/dsgseq/${fas}.count &
  fi

done



#ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq\|fastq" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | cut -d "." -f 1 | sort | uniq | parallel --gnu --will-cite -j $nCores "
for fas in $(ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq\|fastq" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | cut -d "." -f 1 | sort | uniq); do
  echo running htseq...
  
  if [ $pairstate == "paired" ]; then
  
    ! test -d ${outputdir}/results/htseq_exon && mkdir -p ${outputdir}/results/htseq_exon || true

    if ! test -s ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.sam; 
      then
      samtools view -h ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.bam > ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.sam
    fi 

    if ! test -s "${outputdir}/results/htseq_exon/${fas}"* ; then
    #if ! test -d ${outputdir}/results/htseq_exon; then
      echo ${fas} not there
      python ${path}/runner/seqGSEA/count_in_exons.py -p yes ${index}/this_gtf.exon ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.sam ${outputdir} ${fas} ${meta}
    fi
    
    echo DONE HTSeq ${fas}
    

  else 
    echo running for single-end read data
    echo Transcript quantification with Salmon for Sample: ${fas}

    ! test -d ${outputdir}/results/htseq_exon && mkdir -p ${outputdir}/results/htseq_exon || true

    if ! test -s ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.sam; 
      then
      samtools view -h ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.bam > ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.sam

      python ${path}/runner/seqGSEA/count_in_exons.py ${index}/this_gtf.exon ${outputdir}/alignments/${fas}/${fas}Aligned.sortedByCoord.out.sam ${outputdir} ${fas} ${meta}
    fi 

    echo DONE HTSeq ${fas}
    
  fi
done


##### fix the HTSeq count files  ---> remove quotes that generate eroors in seqGSEA
# pushd  ${outputdir}/results/htseq_exon
# for file in $(ls); do sed 's/\"//g' $file > ${file}_clean;  mv ${file}_clean ${file}; done;
# popd


################### DSGSeq ###################
#ls ${readfilesdir} | grep "fasta.gz\|fastq.gz\|fq" | awk -F- '{ print $1 }' | sed 's:.*/::' | sed 's/.\///g' | sed 's:_[^_]*$::' | cut -d "." -f 1 | sort | uniq | parallel --gnu  --will-cite -j $nCores "






if ! test -f ${outputdir}/results/salmon_count.csv;
  then
  echo "Extracting salmon reads..."
  python ${path}/scripts/extract_salmon.py --dir ${outputdir}/salmon_out --outputfile ${outputdir}/results/ --pattern $pattern --type TPM --gtf ${index}/tx2g.gtf --meta $meta
fi

if ! test -f ${outputdir}/results/kal_count.csv;
  then
  echo "Extracting kallisto reads..."
  python ${path}/scripts/extract_kallisto.py --dir ${outputdir}/kallisto_out --outputfile ${outputdir}/results/ --pattern $pattern --type TPM --gtf ${index}/tx2g.gtf --meta $meta
fi

if ! test -f ${outputdir}/results/rsem_count.csv;
  then
  echo "Extracting rsem reads..."
  python ${path}/scripts/extract_rsem.py --dir ${outputdir}/rsem_out --outputfile ${outputdir}/results/ --pattern $pattern --type TPM --gtf $gtf --meta $meta
fi

# conda activate py2
# if ! test -d ${outputdir}/results/cuffdiff_results; then
#   bash ${path}/run_cufflinks.sh --config ${config}
# fi 

# if ! test -d ${outputdir}/results/jseq_${groupl[0]}_${groupl[1]}_; then
#   singularity run --bind ${outputdir}:/MOUNT /nfs/proj/is_benchmark/runner/sing_junctionseq/jcseq.sif 
# fi
#python ${path}/scripts/extract_salmon.py --dir ${readfilesdir}/salmon_out --outputfile ${readfilesdir}/results/salmon_count.csv --pattern sample --type TPM



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
#conda activate nease
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
  java -jar /nfs/scratch/chit/iso-kTSP_v1.0.3.jar ${outputdir}/results/salmon_isoktsp_count.txt -i -n 2 -s 5000 -o ${outputdir}/results/salmon_isoktsp_output_${groupl[0]}_${groupl[1]}.txt --seed 1234 -c ${groupl[0]} ${groupl[1]}
fi

if ! test -f ${outputdir}/results/kal_isoktsp_output_${groupl[0]}_${groupl[1]}.txt;
  then
  java -jar /nfs/scratch/chit/iso-kTSP_v1.0.3.jar ${outputdir}/results/kal_isoktsp_count.txt -i -n 2 -s 5000 -o ${outputdir}/results/kal_isoktsp_output_${groupl[0]}_${groupl[1]}.txt --seed 1234 -c ${groupl[0]} ${groupl[1]}
fi

if ! test -f ${outputdir}/results/rsem_isoktsp_output_${groupl[0]}_${groupl[1]}.txt;
  then
  java -jar /nfs/scratch/chit/iso-kTSP_v1.0.3.jar ${outputdir}/results/rsem_isoktsp_count.txt -i -n 2 -s 5000 -o ${outputdir}/results/rsem_isoktsp_output_${groupl[0]}_${groupl[1]}.txt --seed 1234 -c ${groupl[0]} ${groupl[1]}
fi

echo "Finished running iso-KTSP"
echo "#########################"
echo "Running DEXSeq from R"
echo "#########################"
if ! grep -q "dexseq" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "dexseq" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "dexseq" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/dexseq.R $outputdir $path $meta ${groupl[@]}
fi


R --version
echo "#########################"
echo "Running DTurtle from R"
echo "#########################"
if ! grep -q "dturtle" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "dturtle" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "dturtle" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/dtu.R $outputdir $path $meta $gtf ${groupl[@]}
fi

echo "#########################"
echo "Running iso-KTSP from R"
echo "#########################"
if ! grep -q "iso_ktsp" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "iso_ktsp" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "iso_ktsp" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/iso-ktsp.R $outputdir $path ${groupl[@]}
fi


#Rscript cuffdiff.R


echo "#########################"
echo "Running NBSplice from R"
echo "#########################"
if ! grep -q "nbsplice" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "nbsplice" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "nbsplice" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
  then
  apptainer run --mount type=bind,src=${outputdir},dst=/mnt /nfs/proj/is_benchmark/runner/docker_nbsplice/nbsplice.sif
  # singularity run --bind ${outputdir}:/MOUNT --bind /nfs/proj/is_benchmark/Rlib:/usr/local/lib/R/site-library /nfs/proj/is_benchmark/runner/sing_nbsplice/nbsplice.sif
  #docker run -v ${outputdir}:/MOUNT --rm --name 'nbsplice' nbsplice
  #Rscript ${path}/runner/junctionseq_res.R $outputdir $path $meta ${groupl[@]} 
fi




echo "#########################"
echo "Running Saturn from R"
echo "#########################"
if ! grep -q "saturn" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "saturn" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "saturn" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/saturn.R $outputdir $path $meta ${groupl[@]} 
fi

echo "#########################"
echo "Running DSGseq from R"
echo "#########################"

if ! grep -q "DSGseq" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "DSGseq" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "DSGseq" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; then
  

  ######### change here to fit the group name of meta data ###########
  group1=$(awk -v outfile=${outputdir} 'BEGIN{FS=OFS=" "}{ if ($2 == "T") print outfile"/results/dsgseq/"$1".count" }' ${meta} | grep -v "group") 
  group2=$(awk -v outfile=${outputdir} 'BEGIN{FS=OFS=" "}{ if ($2 == "N") print outfile"/results/dsgseq/"$1".count" }' ${meta} | grep -v "group") 

  pushd /nfs/home/students/chit/DSGseq/DSG-0.1.0/
  Rscript DSGNB.R $(echo ${group1[@]} | wc -w) ${group1[@]} $(echo ${group2[@]} | wc -w) ${group2[@]} ${outputdir}/results/DSGseq_results.diff
  popd

  Rscript ${path}/runner/dsgseq.R $outputdir $path ${groupl[@]}
fi 

echo "#########################"
echo "Running Limma from R"
echo "#########################"
if ! grep -q "LimmaDS" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "LimmaDS" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "LimmaDS" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/limma.R $outputdir $path $meta ${groupl[@]} 
fi

echo "#########################"
echo "Running edgeR from R"
echo "#########################"
if ! grep -q "edgeR" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "edgeR" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "edgeR" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/edgeR.R $outputdir $path $meta ${groupl[@]} 
fi


# echo "#########################"
# echo "Running IUTA from R"
# echo "#########################"

# if ! grep -q "IUTA" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "IUTA" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "IUTA" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
# then
#   cp /nfs/scratch/chit/ref/ensembl_98/iuta_transcriptinfo.gtf ${outputdir}
#   singularity run --bind ${outputdir}:/MOUNT --bind /nfs/proj/is_benchmark/Rlib:/usr/local/lib/R/site-library /nfs/proj/is_benchmark/runner/sing2docker/iuta.sif

# fi


echo "#########################"
echo "Running seqGSEA from R"
echo "#########################"
if ! grep -q "seqGSEA" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "seqGSEA" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "seqGSEA" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; then
  for file in $(ls ${outputdir}/results/htseq_exon); 
  do
    sed 's/"//g' ${outputdir}/results/htseq_exon/${file} > ${outputdir}/results/htseq_exon/temp.txt 
    #mv ${outputdir}/results/htseq_exon/temp.txt ${outputdir}/results/htseq_exon/${file}
  done
  Rscript ${path}/runner/seqGSEA/seqgsea.R $outputdir $path $meta $index ${groupl[@]}
fi


# conda install -c r r r-essentials

echo "#########################"
echo "Running JunctionSeq from R"
echo "#########################"
# if ! grep -q "junctionseq" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "junctionseq" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "junctionseq" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
#   then
#   # singularity run --bind ${outputdir}:/MOUNT --bind /nfs/scratch/chit/is_benchmark/Rlib:/usr/local/lib/R/site-library /nfs/proj/is_benchmark/runner/sing_junctionseq/jcseq.sif 
#   #docker run -v ${outputdir}:/MOUNT --rm --name 'jcseq' jcseq 
#   Rscript ${path}/runner/junctionseq_res.R $outputdir $path $meta ${groupl[@]} 
# fi


echo "#########################"
echo "Running Cuffdiff from R"
echo "#########################"
if ! grep -q "cuffdiff" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "cuffdiff" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "cuffdiff" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  Rscript ${path}/runner/cuffdiff.R $outputdir $path $meta ${groupl[@]} $gtf
fi

echo "#########################"
echo "Running rDiff from R"
echo "#########################"
if ! grep -q "rDiff" ${outputdir}/results/salmon_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "rDiff" ${outputdir}/results/kal_res_gene_${groupl[0]}_${groupl[1]}.txt || ! grep -q "rDiff" ${outputdir}/results/rsem_res_gene_${groupl[0]}_${groupl[1]}.txt; 
then
  singularity run --bind ${outputdir}:/MOUNT --bind /nfs/proj/is_benchmark/Rlib:/usr/local/lib/R/site-library --bind /nfs/scratch/chit/ref/ensembl_98:/ref /nfs/proj/is_benchmark/runner/sing2docker/rdiff.sif 
fi