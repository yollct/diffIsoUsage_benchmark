# install.packages("/nfs/home/students/chit/isoform_0.99.2.tar.gz", type="source", repo=NULL)
suppressMessages(library(tidyverse))
suppressMessages(library(AnnotationDbi))
library(isoform)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)
library(foreach)
library(doParallel)
library(dplyr)
cl <- makeCluster(4)
registerDoParallel(cl)
clusterCall(cl, function() library(magrittr))

path="/nfs/home/students/chit/is_benchmark"
gtf = '/nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.gtf'

# --------------------------------------------------------- 
# cluster by cluster, identify non-overlap exons
# --------------------------------------------------------- 
#### code obtain from isodot isoform package
gtf_df <- gtf_df %>% dplyr::filter(type=="exon") %>% group_by(seqnames, gene_id) %>% mutate(clust=cur_group_id())

gtf_df <- gtf_df %>% group_by(gene_name) %>% mutate(exon_number=1:n()) %>% ungroup

inf <- gtf_df %>% dplyr::filter(type=="exon") %>% dplyr::mutate(clustId=sprintf("chr%s_%s", seqnames, clust)) %>% dplyr::rename(exonId=exon_number, tranId=transcript_ver, tranNm=transcript_name, geneId=gene_id, geneNm=gene_name) %>%  dplyr::select(seqnames, start, end, clustId, exonId, tranId, tranNm, geneId, geneNm, source, type, strand, score, phase) %>% arrange(seqnames, start) 


infN  = matrix("", nrow=1000000, ncol=14)

ucIds = unique(inf$clustId)
length(ucIds)

kk   = 0
idx1 = idx2 = 1

pasteUniqu = function(v){
  it = unlist(strsplit(v, ":", fixed=TRUE))
  paste(sort(unique(it)),collapse=":")
}



for(i in 1:length(ucIds)){ 
  kk = kk + 1
  if(kk %% 1000 == 0){
    cat(kk, date(), "\n")
  }
  
  ## inf1 is the annotation for one cluster ID
  inf1 = inf %>% dplyr::filter(clustId==ucIds[i]) %>% dplyr::arrange(start, end)
  nn   = nrow(inf1)
  
  
  ## if this cluster only has one exon
  if(nn == 1) { 
    inf1$exonId = 1
    infN[idx1,] = as.matrix(inf1)
    idx1 = idx1 + 1
    next 
  }
  
  ## calcualte gaps
  gaps1   = inf1$start[-1] - inf1$end[-nn]  
  w2check = which(gaps1 <= 0)
    
  ## if there is no overlapping exons
  if(length(w2check) == 0){
    inf1$exonId = 1:nn
    idx2 = idx1 + nrow(inf1) - 1
    infN[idx1:idx2,] = as.matrix(inf1)
    idx1 = idx2 + 1
    next 
  }
    
  while(length(w2check) > 0){
    
    pos1 = w2check[1]
    pos2 = pos1+1
    
    ## well, we start with unique exons, but after a while...
    if(inf1$end[pos1] == inf1$end[pos2] && inf1$start[pos1] == inf1$start[pos2]){
      ## if two exons share both the start and the end poistions
      inf1$tranId[pos1] = pasteUniqu(inf1$tranId[pos1:pos2])
      inf1$tranNm[pos1] = pasteUniqu(inf1$tranNm[pos1:pos2])
      inf1$geneId[pos1] = pasteUniqu(inf1$geneId[pos1:pos2])
      inf1$geneNm[pos1] = pasteUniqu(inf1$geneNm[pos1:pos2])
      inf1              = inf1[-pos2,]
      
    }else if(inf1$end[pos1] == inf1$end[pos2]){
      ## if two exons share the end poistions

      ## if the start position difference is no greater than 2, 
      ## combine them
      if(abs(inf1$start[pos1] - inf1$start[pos2]) < 3){
        inf1$tranId[pos1] = pasteUniqu(inf1$tranId[pos1:pos2])
        inf1$tranNm[pos1] = pasteUniqu(inf1$tranNm[pos1:pos2])
        inf1$geneId[pos1] = pasteUniqu(inf1$geneId[pos1:pos2])
        inf1$geneNm[pos1] = pasteUniqu(inf1$geneNm[pos1:pos2])
        inf1              = inf1[-pos2,]
      }else{
        inf1$end[pos1]    = inf1$start[pos2] - 1
        inf1$tranId[pos2] = pasteUniqu(inf1$tranId[pos1:pos2])
        inf1$tranNm[pos2] = pasteUniqu(inf1$tranNm[pos1:pos2])
        inf1$geneId[pos2] = pasteUniqu(inf1$geneId[pos1:pos2])
        inf1$geneNm[pos2] = pasteUniqu(inf1$geneNm[pos1:pos2])
      }
    }else if(inf1$start[pos1] == inf1$start[pos2]){
      ## if two exons share the start poistions
      
      ## if the end position difference is no greater than 2, 
      ## combine them
      if(abs(inf1$end[pos1] - inf1$end[pos2]) < 3){
        inf1$tranId[pos1] = pasteUniqu(inf1$tranId[pos1:pos2])
        inf1$tranNm[pos1] = pasteUniqu(inf1$tranNm[pos1:pos2])
        inf1$geneId[pos1] = pasteUniqu(inf1$geneId[pos1:pos2])
        inf1$geneNm[pos1] = pasteUniqu(inf1$geneNm[pos1:pos2])
        inf1              = inf1[-pos2,]
      }else{
        inf1$start[pos2]  = inf1$end[pos1] + 1
        inf1$tranId[pos1] = pasteUniqu(inf1$tranId[pos1:pos2])
        inf1$tranNm[pos1] = pasteUniqu(inf1$tranNm[pos1:pos2])
        inf1$geneId[pos1] = pasteUniqu(inf1$geneId[pos1:pos2])
        inf1$geneNm[pos1] = pasteUniqu(inf1$geneNm[pos1:pos2])
      }
    }else if(inf1$end[pos2] < inf1$end[pos1]){
      ## if the 2nd exon is within the first exon
      newExon           = inf1[pos1,,drop=FALSE]
      newExon$start     = inf1$end[pos2]+1
      newExon$end       = inf1$end[pos1]
                
      inf1$end[pos1]    = inf1$start[pos2] - 1
      
      inf1$tranId[pos2] = pasteUniqu(inf1$tranId[pos1:pos2])
      inf1$tranNm[pos2] = pasteUniqu(inf1$tranNm[pos1:pos2])
      inf1$geneId[pos2] = pasteUniqu(inf1$geneId[pos1:pos2])
      inf1$geneNm[pos2] = pasteUniqu(inf1$geneNm[pos1:pos2])

      nn = nrow(inf1)
      ## new exon should be after pos2, but inf1 will be re-ordered anyway
      inf1 = rbind(inf1[1:pos1,,drop=FALSE], newExon, inf1[pos2:nn,,drop=FALSE])
      
    }else{
      
      newExon          = inf1[pos1,,drop=FALSE]
      newExon$start    = inf1$start[pos2]
      newExon$end      = inf1$end[pos1]
      newExon$tranId   = pasteUniqu(inf1$tranId[pos1:pos2])
      newExon$tranNm   = pasteUniqu(inf1$tranNm[pos1:pos2])
      newExon$geneId   = pasteUniqu(inf1$geneId[pos1:pos2])
      newExon$geneNm   = pasteUniqu(inf1$geneNm[pos1:pos2])

      inf1$end[pos1]   = newExon$start - 1
      inf1$start[pos2] = newExon$end + 1
      
      nn = nrow(inf1)
      inf1 = rbind(inf1[1:pos1,,drop=FALSE], newExon, inf1[pos2:nn,,drop=FALSE])
      
    }
    

    nn = nrow(inf1)
    if(nn==1){ break } 
    
    if(any(inf1$end < inf1$start)){
      stop("something is wrong...\n")
    }
    
    inf1 = inf1[order(inf1$start, inf1$end),]
    
    gaps1 = inf1$start[-1] - inf1$end[-nn]  
    w2check = which(gaps1 <= 0)
  }
  
  inf1$exonId = 1:nn

  idx2 = idx1 + nrow(inf1) - 1
  infN[idx1:idx2,] = as.matrix(inf1)
  idx1 = idx2 + 1

}

infN[(idx1-1):idx1,]
infN = infN[1:(idx1-1),]

infN = as.data.frame(infN, stringsAsFactors=FALSE)
write.table()


colnames(infN) = c(names(inf))

infN$exonId = as.numeric(infN$exonId)

infN <- readRDS("/nfs/data/covid_hscell_4tp/ensembl_106/infN.rds")
write.table(infN, col.names = FALSE, append = FALSE, file = "/nfs/data/covid_hscell_4tp/ensembl_106/hg38_nonoverlapexon.bed", quote = FALSE, sep = "\t", row.names = FALSE)

geneId2use  = paste("gene_id \"", infN$geneId, "\";", sep="")
tranId2use  = paste("transcript_id \"", infN$tranId, "\";", sep="")
geneNm2use  = paste("gene_name \"", infN$geneNm, "\";", sep="")
tranNm2use  = paste("transcript_name \"", infN$tranNm, "\";", sep="")
exonId2use  = paste("exon_id \"", infN$exonId, "\";", sep="")
clustId2use = paste("clustId \"", infN$clustId, "\";", sep="")

infN$anno   = paste(geneId2use, tranId2use, geneNm2use, 
                    tranNm2use, exonId2use, clustId2use, sep=" ")

dim(infN)
infN[1:2,]
infN[infN==""] = NA
# --------------------------------------------------------- 
# double check the gaps between exons
# --------------------------------------------------------- 

nn = nrow(infN)
infN$start = as.numeric(infN$start)
infN$end   = as.numeric(infN$end)

gaps1 = infN$start[-1] - infN$end[-nn]
gaps2 = infN$start[-1] - infN$start[-nn]

w2ck = which(infN$clustId[-1] == infN$clustId[-nn])

gaps1 = gaps1[w2ck]
gaps2 = gaps2[w2ck]

summary(gaps1)
summary(gaps2)



# png("ensembl_gaps_nonoverlap_exons.png", width=12, height=8, units="in", res=200)
# par(mar=c(5,4,2,1), mfrow=c(2,1))

# hist(gaps1[abs(gaps1)<200], breaks=40, ylab="Frequency", main="",
# xlab="Start of the (i+1)-th exon - End of the i-th exon")
# hist(gaps2[abs(gaps2)<200], breaks=40, ylab="Frequency", main="",
# xlab="Start of the (i+1)-th exon - Start of the i-th exon")

# dev.off()

# len = infN$end - infN$start
# summary(len)
# table(len==0)

# png("ensembl_len_nonoverlap_exons.png", width=5, height=4, units="in", res=200)
# par(mar=c(5,4,2,1))
# hist(log10(len), ylab="Frequency", main="",
# xlab="Exon length, log10(bp)")

# dev.off()

# ---------------------------------------------------------
# write out results
# ---------------------------------------------------------



write.table(infN[, c(1:8,ncol(infN))], col.names = FALSE, append = FALSE, 
  file = "/nfs/data/covid_hscell_4tp/ensembl_106/hg38_exon_isodot.bed", 
  quote = FALSE, sep = "\t", row.names = FALSE)
