library(IUTA)
library(dplyr)

args <- commandArgs(trailingOnly=T)
outdir <- args[1]
con1 <- args[2]
con2 <- args[3]

metadata <- read.csv("/MOUNT/meta.txt", sep="\t")
group1 <- metadata %>% dplyr::filter(group==con1)
group2 <- metadata %>% dplyr::filter(group==con2)

bamdir <- "/MOUNT/alignments/"
bam.list.1 <- lapply(group1$sample_id, function(x){paste0(bamdir,"/",x,"/",x, "Aligned.sortedByCoord.out.bam")}) %>% unlist
bam.list.2 <- lapply(group2$sample_id, function(x){paste0(bamdir,"/",x,"/",x, "Aligned.sortedByCoord.out.bam")}) %>% unlist
print("check")
bam.list.1
print(length(bam.list.1))
starttime <- Sys.time()
##
IUTA(
    bam.list.1=bam.list.1,
    bam.list.2=bam.list.2,
    transcript.info="/MOUNT/iuta_transcriptinfo.gtf",
    rep.info.1= rep(1, length(bam.list.1)),
    rep.info.2= rep(1,length(bam.list.2)),
    output.dir="/MOUNT/results/",
    output.na=FALSE,
    FLD="normal"
)
endtime <- Sys.time()

print("IUTA took ")
print(starttime-endtime)
