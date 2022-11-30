suppressMessages(library(tidyverse))
suppressMessages(library(AnnotationDbi))
#install.packages("/nfs/home/students/chit/isoform_0.99.2.tar.gz", repos=NULL, type="source")
library(isoform)
library(rtracklayer)
library(GenomicRanges)
library(foreach)
library(doParallel)
library(dplyr)
cl <- makeCluster(10)
registerDoParallel(cl)
clusterCall(cl, function() library(magrittr))
path <- "/nfs/home/students/chit/is_benchmark"

meta1 = read.csv(paste0(path,"/simulated_reads/sim_rep_info.txt"), sep="\t")
colnames(meta1) <- c("sample_id", "group", "lib_sizes")
meta1

idx1 = 1
idx2 = 0
tags = as.character(meta1$sample_id)
Routs = paste(path,"/alignments/",tags,"/",tags,"_isodot_output_knowniso.RData", sep="")
fragSizeFiles = paste(path,"/alignments/",tags,"/",tags,"_insertSize_dist.txt", sep="")

countdepth = paste(path,"/alignments/",tags,"/",tags,"_counts_read_isodot.txt", sep="")
readDep <- do.call("c", lapply(countdepth, function(x){
    one = readLines(countdepth[1])
    allcounts <- do.call("sum", lapply(1:length(one), function(x){
        tmp <- strsplit(one[x], " ")
        tmp <- tmp[[1]][tmp[[1]]!=""]
        as.numeric(tmp[1])
    }))
    allcounts
}))



# xData is a case/control indicator of the two samples
xData = meta1$group
# -------------------------------------------------------------------------
# select the transcript clusters to be used
# -------------------------------------------------------------------------
allct <- c()
for (i in 1:length(Routs)){
    load(Routs[i])
    ct1 = sapply(geneMod, function(x){ sum(x$y) } )
    allct[[i]]<-names(ct1)
}

g2test = Reduce(intersect, allct)
length(g2test)
# -------------------------------------------------------------------------
# the cases where there are reads from one sample but not the other sample
# are dropped from testing, though they may be very intersting results.
# -------------------------------------------------------------------------
g2drop = setdiff(Reduce(union, allct), g2test)
c2drop = data.frame(tcluster=g2drop, ct1=ct1[match(g2drop,names(ct1))])

# # -------------------------------------------------------------------------
# # test differential expression and differential isoform usage
# # -------------------------------------------------------------------------
# outputFileName = "set1_vs_set2_knownIsoforms"
# outputFileName = sprintf("%s_%d_%d.txt", outputFileName, idx1, idx2)
# isoDu(tags, Routs, xData, outputFileName, fragSizeFiles, g2test=g2test,
# readLen=76, method="bootstrap", lmax=500, idxs=idx1:idx2, maxTime=43200)
gL <- c()
for (i in 1:length(Routs)){
    load(Routs[i])
    gL[[i]] <- geneMod
}

fragsize <- readFragmentSizes(fragSizeFiles)
# -------------------------------------------------------------------------
# test differential isoform usage only
# -------------------------------------------------------------------------
outputFileName = paste0(path,"/results/isodot_dtu")
outputFileName = sprintf("%s_%d_%d.txt", outputFileName, idx1, idx2)
isoDu(tags, gL, xData, readDep, outputFileName, fragsize, g2test=g2test,
readLen=76, method="bootstrap", lmax=500, duOnly=TRUE,
maxTime=43200)

