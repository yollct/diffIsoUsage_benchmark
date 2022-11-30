library(rtracklayer)
library(rtracklayer)

path <- "/nfs/home/students/chit/is_benchmark"
gtf <-rtracklayer::import(paste0(path,"/simulation/splicing_variants.gtf"))

tr_gtf <- as.data.frame(gtf) %>% dplyr::filter(type=="transcript")
tr_gtf %>% filter(template==TRUE)
rtracklayer::export(tr_gtf, paste0(path,"/simulation/splicing_variants_transcript.gtf"), format="gtf")
