library(NBSplice)
library(tximport)
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
meta <- args[2]
con1 <- args[3]
con2 <- args[4]
con <- c(con1, con2)
print(con)

meta1 = read.csv(meta, sep="\t")
row.names(meta1) <- meta1$sample_id

meta1 <- meta1 %>% dplyr::filter(group %in% con)
meta1$sample_id <- lapply(meta1$sample_id, function(x){gsub(".sra", "", x)}) %>% unlist

resgene <- read.csv(paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), sep="\t")
restx <- read.csv(paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), sep="\t")
# sresgene <- read.csv(paste0(outdir, sprintf("/results/stager_salmon_res_gene_%s_%s.txt", con1, con2)), sep="\t")


if (!any(grepl("nbsplice", colnames(resgene)))) { 
    genename  <- read.csv(paste0(outdir, "/results/salmon_count.csv"))
    files <- Sys.glob(paste0(outdir, "/salmon_out/*/quant.sf"))
    names(files) <- gsub(".*/","",gsub("/quant.sf","",files))
    txi <- tximport(files, type="salmon", txOut=TRUE)
    isoCounts <- txi$counts
    row.names(isoCounts) <- lapply(row.names(isoCounts), function(x){strsplit(x,"[.]")[[1]][1]}) %>% unlist
    print(genename %>% head)
    geneIso <- data.frame(gene_id=genename$gene_id, isoform_id=genename$feature_id)
    geneIso$isoform_id <- lapply(geneIso$isoform_id, function(x){strsplit(x,"[.]")[[1]][1]}) %>% unlist
    

    designMatrix <- meta1 %>% dplyr::select(sample_id, group) %>% dplyr::rename(sample=sample_id, condition=group)
    row.names(designMatrix) <- designMatrix$sample
    designMatrix$condition <- factor(designMatrix$condition, levels=con)
    
    print(row.names(designMatrix) %in% colnames(isoCounts))
    print(dim(isoCounts))
    isoCounts <- as.data.frame(isoCounts)
    isoCounts <- isoCounts %>% dplyr::select(row.names(designMatrix))

    colName <- "condition"
    myIsoDataSet<-IsoDataSet(isoCounts, designMatrix, colName, geneIso)

    myIsoDataSet<-buildLowExpIdx(myIsoDataSet, colName, ratioThres = 0.01, 
                            countThres = 1)

    myDSResults<-NBTest(myIsoDataSet, colName, test="F")

    res <- results(myDSResults) 
    nbgene <- res %>% dplyr::select(gene, geneFDR) %>% group_by(gene) %>% summarise(nbsplice=min(geneFDR)) %>% dplyr::rename(feature_id=gene)
    resgene <- full_join(resgene, nbgene, by="feature_id")

    nbtx <- res %>% dplyr::select(iso, FDR) %>% dplyr::rename(feature_id=iso)
    restx <- inner_join(restx, nbtx, by="feature_id")

    write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
    
    ##stageR 
    # pConfirmation <- matrix(res$pval)
    # dimnames(pConfirmation) <- list(c(res$iso), c("transcript"))

    # qval <- res$genePval
    # qval[qval==0] <- min(qval[qval !=0]) * 0.1
    # pScreen <- qval
    # tx2gene <- res[,c("iso", "gene")]
    # stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene) 
    # stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=T) 
    # suppressWarnings({ padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE) }) 

    # stagegene <- padj %>% dplyr::select(geneID, gene) %>% unique()
    # stagegene <- stagegene %>% dplyr::rename(nbsplice_stageR=gene, feature_id=geneID)
    # sresgene <- full_join(stagegene, sresgene, by=c("feature_id"="feature_id"))

    # stagetx <- padj %>% dplyr::select(txID, transcript) %>% unique()
    # stagetx <- stagetx %>% dplyr::rename(nbsplice_stageR=transcript, feature_id=txID)
    # stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    # srestx <- full_join(stagetx, srestx, by=c("feature_id"="feature_id"))

    # write.table(sresgene, paste0(outdir, sprintf("/results/stager_salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    # write.table(srestx, paste0(outdir, sprintf("/results/stager_salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}

rm(genename)
rm(salmoncnt)
rm(cnt)
rm(d)
rm(resgene)
rm(restx)
rm(sresgene)
rm(srestx)

resgene <- read.csv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), sep="\t")
restx <- read.csv(paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), sep="\t")

# sresgene <- read.csv(paste0(outdir, sprintf("/results/stager_kal_res_gene_%s_%s.txt", con1, con2)), sep="\t")



if (!any(grepl("nbsplice", colnames(resgene)))) {
    print("kallisto count nbsplice")
    genename  <- read.csv(paste0(outdir, "/results/kal_count.csv"), sep="\t")
    files <- Sys.glob(paste0(outdir, "/kallisto_out/*/abundance.h5"))
    names(files) <- gsub(".*/","",gsub("/abundance.h5","",files))
    txi <- tximport(files, type="kallisto", txOut=TRUE)
    isoCounts <- txi$counts
    row.names(isoCounts) <- lapply(row.names(isoCounts), function(x){strsplit(x,"[.]")[[1]][1]}) %>% unlist

    geneIso <- data.frame(gene_id=genename$gene_id, isoform_id=genename$feature_id)
    geneIso$isoform_id <- lapply(geneIso$isoform_id, function(x){strsplit(x,"[.]")[[1]][1]}) %>% unlist

    
    designMatrix <- meta1 %>% dplyr::select(sample_id, group) %>% dplyr::rename(sample=sample_id, condition=group)
    row.names(designMatrix) <- designMatrix$sample
    designMatrix$condition <- factor(designMatrix$condition, levels=con)

    isoCounts <- as.data.frame(isoCounts)
    isoCounts <- isoCounts %>% dplyr::select(row.names(designMatrix))

    colName <- "condition"
    myIsoDataSet<-IsoDataSet(isoCounts, designMatrix, colName, geneIso)

    myIsoDataSet<-buildLowExpIdx(myIsoDataSet, colName, ratioThres = 0.01, 
                            countThres = 1)

    myDSResults<-NBTest(myIsoDataSet, colName, test="F")

    res <- results(myDSResults) 
    nbgene <- res %>% dplyr::select(gene, geneFDR) %>% group_by(gene) %>% summarise(nbsplice=min(geneFDR)) %>% dplyr::rename("feature_id"="gene")
    print(nbgene %>% head)
    resgene <- full_join(resgene, nbgene, by="feature_id")

    nbtx <- res %>% dplyr::select(iso, FDR) %>% dplyr::rename(feature_id=iso)
    restx <- inner_join(restx, nbtx, by="feature_id")

    write.table(resgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
    
    ##stageR 
#     pConfirmation <- matrix(res$pval)
#     dimnames(pConfirmation) <- list(c(res$iso), c("transcript"))

#     qval <- res$genePval
#     qval[qval==0] <- min(qval[qval !=0]) * 0.1
#     pScreen <- qval
#     tx2gene <- res[,c("iso", "gene")]
#     stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene) 
#     stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=T) 
#     suppressWarnings({ padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE) }) 

#     stagegene <- padj %>% dplyr::select(geneID, gene) %>% unique()
#     stagegene <- stagegene %>% dplyr::rename(nbsplice_stageR=gene, feature_id=geneID)
#     sresgene <- full_join(stagegene, sresgene, by=c("feature_id"="feature_id"))

#     stagetx <- padj %>% dplyr::select(txID, transcript) %>% unique()
#     stagetx <- stagetx %>% dplyr::rename(nbsplice_stageR=transcript, feature_id=txID)
#     stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
#     srestx <- full_join(stagetx, srestx, by=c("feature_id"="feature_id"))

#     write.table(sresgene, paste0(outdir, sprintf("/results/stager_kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
#     write.table(srestx, paste0(outdir, sprintf("/results/stager_kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}

rm(genename)
rm(salmoncnt)
rm(cnt)
rm(d)
rm(resgene)
rm(restx)
rm(sresgene)
rm(srestx)

resgene <- read.csv(paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), sep="\t")

restx <- read.csv(paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)), sep="\t")
# sresgene <- read.csv(paste0(outdir, sprintf("/results/stager_rsem_res_gene_%s_%s.txt", con1, con2)), sep="\t")


if (!any(grepl("nbsplice", colnames(resgene)))) {
    print("running rsem count nbsplice")
    files <- Sys.glob(paste0(outdir, "/rsem_out/*/*.isoforms.results"))
    names(files) <- gsub(".*/","",gsub("/*.isoforms.results","",files))
    genename  <- read.csv(files[1], sep="\t")
    txi <- tximport(files, type="rsem", txOut=TRUE)
    isoCounts <- txi$counts
    row.names(isoCounts) <- lapply(row.names(isoCounts), function(x){strsplit(x,"[.]")[[1]][1]}) %>% unlist


    geneIso <- data.frame(gene_id=genename$gene_id, isoform_id=genename$transcript_id)
    

    designMatrix <- meta1 %>% dplyr::select(sample_id, group) %>% dplyr::rename(sample=sample_id, condition=group)
    row.names(designMatrix) <- designMatrix$sample
    designMatrix$condition <- factor(designMatrix$condition, levels=con)

    isoCounts <- as.data.frame(isoCounts)
    isoCounts <- isoCounts %>% dplyr::select(row.names(designMatrix))

    colName <- "condition"
    myIsoDataSet<-IsoDataSet(isoCounts, designMatrix, colName, geneIso)

    myIsoDataSet<-buildLowExpIdx(myIsoDataSet, colName, ratioThres = 0.01, 
                            countThres = 1)

    myDSResults<-NBTest(myIsoDataSet, colName, test="F")

    res <- results(myDSResults) 
    nbgene <- res %>% dplyr::select(gene, geneFDR) %>% group_by(gene) %>% summarise(nbsplice=min(geneFDR)) %>% dplyr::rename(feature_id=gene)

    resgene <- full_join(resgene, nbgene, by="feature_id")

    nbtx <- res %>% dplyr::select(iso, FDR) %>% dplyr::rename(feature_id=iso)
    restx <- inner_join(restx, nbtx, by="feature_id")

    write.table(resgene, paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
    
    ##stageR 
    # pConfirmation <- matrix(res$pval)
    # dimnames(pConfirmation) <- list(c(res$iso), c("transcript"))

    # qval <- res$genePval
    # qval[qval==0] <- min(qval[qval !=0]) * 0.1
    # pScreen <- qval
    # tx2gene <- res[,c("iso", "gene")]
    # stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene) 
    # stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=T) 
    # suppressWarnings({ padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE) }) 

    # stagegene <- padj %>% dplyr::select(geneID, gene) %>% unique()
    # stagegene <- stagegene %>% dplyr::rename(nbsplice_stageR=gene, feature_id=geneID)
    # sresgene <- full_join(stagegene, sresgene, by=c("feature_id"="feature_id"))

    # stagetx <- padj %>% dplyr::select(txID, transcript) %>% unique()
    # stagetx <- stagetx %>% dplyr::rename(nbsplice_stageR=transcript, feature_id=txID)
    # stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    # srestx <- full_join(stagetx, srestx, by=c("feature_id"="feature_id"))

    # write.table(sresgene, paste0(outdir, sprintf("/results/stager_rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    # write.table(srestx, paste0(outdir, sprintf("/results/stager_rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}