library(tidyverse)
library(AnnotationDbi)

path <- "/nfs/home/students/chit/is_benchmark"
source(paste0(path,"/runner/groundtruth.R"))

txcount = read.csv(paste0(path,"/simulated_reads/isoktsp_count.txt"), sep="\t")
txcount$X <- row.names(txcount)

txcount$gene <- lapply(txcount$X, function(x){strsplit(x, ",")[[1]][1]}) %>% unlist
txcount$isoform <- lapply(txcount$X, function(x){strsplit(x, ",")[[1]][2]}) %>% unlist
txcount <- txcount %>% dplyr::select(-X)

txcount <- txcount %>% dplyr::filter(gene!="?") %>% arrange(gene)
ntx <- txcount %>% group_by(gene) %>% summarise(countx=n())
genesum <- txcount %>% dplyr::select(-isoform) %>% group_by(gene) %>% summarise(across(colnames(txcount)[grepl("N", colnames(txcount))],max), across(colnames(txcount)[grepl("T", colnames(txcount))],max), .groups="keep")

genesum_mat_txcount <- genesum[rep(seq_len(dim(genesum)[1]), ntx$countx),]

res <- as.matrix(txcount[colnames(txcount)[grepl("sample", colnames(txcount))]])/as.matrix(genesum_mat_txcount[colnames(genesum_mat_txcount)[grepl("sample",colnames(genesum_mat_txcount))]])

res[is.nan(res)]<-0
res[is.infinite(res)]<-0
eventim <- data.frame(eventimportance = rowMeans(res))
eventim$isoform <- txcount$isoform
eventim$gene <- txcount$gene



filter_eventim <- function(cutoff=0.1, type="gene"){
    if (type=="gene"){
        return(eventim %>% dplyr::filter(eventimportance>cutoff) %>% dplyr::select(gene))
    } else {
        return(eventim %>% dplyr::filter(eventimportance>cutoff) %>% dplyr::select(isoform))
    }
}

eventim %>% filter(gene=="ENSG00000225618")
eventim %>% head

txcount %>% filter(gene=="ENSG00000225618")

