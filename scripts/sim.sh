#!/bin/bash

rm -rf ../simulation
Rscript asim.R -o /nfs/home/students/chit/is_benchmark/simulation -g 1000
Rscript gtf2bed.R

cd ../simulation
convert2bed -i gtf --attribute-key="transcript_id" < splicing_variants_transcript.gtf > splicing_variants.bed


cd ../annotations
/nfs/home/students/chit/bedtools2/bin/fastaFromBed -fi Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed /nfs/home/students/chit/is_benchmark/simulation/splicing_variants.bed -nameOnly -split > splicing_variants.fa 