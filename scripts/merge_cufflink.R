library(dplyr)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)

path <- "/nfs/home/students/chit/is_benchmark/"
samples <- list.files(paste0(path,"alignments/"))

###RSEM
collist <-c()
for (x in samples){
  tmp <- read.csv(paste0(path, "alignments/", x, "/quant.isoforms.results"), sep="\t")
  
  tmplist <- tmp %>% dplyr::select(transcript_id, expected_count)
  colnames(tmplist) <- c("transcript_id", x)
  if (length(collist) == 0){
    collist <- tmplist
  } else {
    collist <- inner_join(collist, tmplist, by="transcript_id")
  }
}

write.table(collist, paste0(path, "simulation/rsem_isocounts.csv"), row.names = T, sep="\t")

##Cufflink
cuffcom <- read.csv(paste0(path, "scripts/merged_asm/cuffcmp.merged.gtf.refmap"), sep="\t")
gtf <- rtracklayer::import(paste0(path, "scripts/merged_asm/merged.gtf")) 
gtf <- as.data.frame(gtf) #%>% dplyr::filter(type=="transcript")
gtf$refmap_id <- paste0(gtf$gene_id,"|",gtf$transcript_id)
mapper <- unique(gtf %>% dplyr::select(refmap_id, oId))
refmap <- inner_join(cuffcom, mapper, by=c("cuff_id_list"="refmap_id"))
mapper_func <- function(x){
  if (grepl("CUFF", x)){
    return(refmap[refmap$oId==x,]$ref_id[1])
  } else {
    return(x)
  }
}

collist <-c()
for (x in samples){
  tmp <- read.csv(paste0(path, "alignments/", x, "/isoforms.fpkm_tracking"), sep="\t")

  tmplist <- tmp %>% dplyr::select(tracking_id, FPKM)
  colnames(tmplist) <- c("transcript_id", x)
  if (length(collist) == 0){
    collist <- tmplist
  } else {
    collist <- inner_join(collist, tmplist, by="transcript_id")
  }
  gc()
}

write.table(collist, paste0(path, "simulated_reads/cufflink_isocounts.csv"), row.names = T, sep="\t")
