# install.packages("/nfs/home/students/chit/isoform_0.99.2.tar.gz", type="source", repo=NULL)
suppressMessages(library(tidyverse))
suppressMessages(library(AnnotationDbi))
library(isoform)
library(rtracklayer)
library(GenomicRanges)
install.packages("/nfs/home/students/chit/plyranges_1.16.0.tar.gz", repos=NULL, type="source")
library(plyranges)
library(foreach)
library(doParallel)
library(dplyr)
cl <- makeCluster(10)
registerDoParallel(cl)
clusterCall(cl, function() library(magrittr))

args = commandArgs(trailingOnly=TRUE)
path="/nfs/home/students/chit/is_benchmark"
gtf = '/nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.gtf'

## things we need
bedfile = "/nfs/data/covid_hscell_4tp/ensembl_106/nonoverlapexon.txt"
sample = args[1]
#sample="sample_02"
bamfile = sprintf("./alignments/%s/%sAligned.sortedByCoord.out.bam", sample, sample)
ctfile = sprintf("./alignments/%s/%s_counts_read_isodot.txt", sample, sample)
fragsize = sprintf("./alignments/%s/%s_insertSize_dist.txt", sample, sample)
output = sprintf("./alignments/%s/%s_isodot_output_knowniso.RData", sample, sample)


print(sprintf("Count reads for %s", sample))
# count the number of reads
countReads(bamfile, bedfile, ctfile)

# isoAll <- knownIsoforms("/nfs/data/covid_hscell_4tp/ensembl_106/hg38_exon_isodot.bed")
# save(isoAll, file="/nfs/data/covid_hscell_4tp/ensembl_106/hg38_knownisoforms.RData")
iso <- "/nfs/data/covid_hscell_4tp/ensembl_106/hg38_knownisoforms.RData"

isoDetector(ctfile, bedfile, fragsize, 76, output, knownIsoforms=iso)





