suppressMessages(library(DEXSeq))
suppressMessages(library(tximport))
suppressMessages(library(DRIMSeq))
suppressMessages(library(tidyverse))
suppressMessages(library(stageR))
suppressMessages(library(rhdf5))
suppressMessages(library(DTUrtle))
biocpar <- BiocParallel::MulticoreParam(12)

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
path <- args[2]
meta <- args[3]
con1 <- args[4]
con2 <- args[5]
con <- c(con1, con2)
print(con)


print("in dexseq")
resgene <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)))

# sresgene <- read_tsv(paste0(outdir, sprintf("/results/stager_salmon_res_gene_%s_%s.txt", con1, con2)))
# srestx <- read_tsv(paste0(outdir, sprintf("/results/stager_salmon_res_tx_%s_%s.txt", con1, con2)))

meta1 = read.csv(meta, sep="\t")
row.names(meta1) <- meta1$sample_id
meta1$sample_id <- lapply(meta1$sample_id, function(x){gsub(".sra", "",x)}) %>% unlist

if (!any(grepl("dexseq", colnames(resgene)))) { 
    genename  <- read.csv(paste0(outdir, "/results/salmon_count.csv"))
    files <- Sys.glob(paste0(outdir, "/salmon_out/*/quant.sf"))
    names(files) <- gsub(".*/","",gsub("/quant.sf","",files))
    txi <- tximport(files, type="salmon", txOut=TRUE)
    salmontpm <- txi$counts

    salmoncnt <- data.frame(gene_id=genename$gene_id, feature_id=genename$feature_id, txi$counts)

    # load(paste0(outdir,"/simulated_reads/sim_counts_matrix.rda"))
    # row.names(counts_matrix) <- lapply(row.names(counts_matrix), function(x){strsplit(x," ")[[1]][1]})

    # row.names(counts_matrix)[!row.names(counts_matrix) %in% row.names(salmoncnt)]
    # counts_matrix[row.names(counts_matrix)=="ENST00000366839.8",]

    salmontpm <- salmontpm[rowSums(salmontpm)>10,]
    salmontpm <- salmontpm[salmoncnt$gene_id != "?"]
    salmoncnt <- salmoncnt %>% dplyr::filter(gene_id!="?")
    colnames(salmoncnt) <- lapply(colnames(salmoncnt), function(x){gsub(".fq", "", x)}) %>% unlist

    meta1 <- meta1 %>% dplyr::filter(group %in% con)

    ## 
    salmoncnt <- salmoncnt %>% dplyr::select("gene_id", "feature_id", meta1$sample_id)


    d <- dmDSdata(counts=salmoncnt, samples=meta1)

    dxd <- DEXSeqDataSet(countData=round(as.matrix(DRIMSeq::counts(d)[,-c(1:2)])), 
                        sampleData=DRIMSeq::samples(d),
                        design=~sample+exon+group:exon,
                        featureID=counts(d)$feature_id,
                        groupID=counts(d)$gene_id)

    dxd <- estimateSizeFactors(dxd)
    dxd <- estimateDispersions(dxd, BPPARAM = biocpar, quiet=T)
    dxd <- testForDEU(dxd, BPPARAM = biocpar, reducedModel = ~sample+exon)
    dxd = estimateExonFoldChanges( dxd, fitExpToVar="group", BPPARAM = biocpar)
    dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
    qval <- perGeneQValue(dxr)
    dxr <- as.data.frame(dxr)
    qval[qval==0] <- min(qval[qval !=0]) * 0.1
    dxr.g <- data.frame(gene=names(qval), qval)
    # dxr$pvalue[dxr$pvalue==0] <- min(dxr$pvalue[dxr$pvalue!=0]) * 0.1


    columns <- c("featureID","groupID","pvalue", "padj") 
    dxr <- as.data.frame(dxr[,columns])

    pConfirmation <- matrix(dxr$pvalue,ncol=1)
    dimnames(pConfirmation) <- list(c(dxr$featureID),c("transcript"))
    pScreen <- qval
    tx2gene <- counts(d)[,c("feature_id", "gene_id")]


    stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene) 
    stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=T) 
    suppressWarnings({ dex.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE) }) 

    
    #resgene <- inner_join(resgene, dxr %>% dplyr::select(groupID, pvalue), by=c("feature_id"="groupID"))
    dxr.g <- dxr.g %>% unique()
    resgene <- full_join(resgene, dxr.g, by=c("feature_id"="gene"))
    resgene <- resgene %>% dplyr::rename(`dexseq`=qval)

    # stagegene <- dex.padj %>% dplyr::select(geneID, gene) %>% unique()
    # stagegene <- stagegene %>% dplyr::rename(dexseq_stageR=gene, feature_id=geneID)
    # sresgene <- full_join(stagegene, sresgene, by=c("feature_id"="feature_id"))
    # #resgene <- resgene %>% dplyr::rename(`dexseq_stageR`=gene, feature_id=geneID)


    dxr.tx <- dxr %>% dplyr::select(featureID, padj) %>% dplyr::rename(dexseq=padj, feature_id=featureID)
    dxr.tx$feature_id <- lapply(dxr.tx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, dxr.tx, by=c("feature_id"="feature_id"))

    # stagetx <- dex.padj %>% dplyr::select(txID, transcript) %>% unique()
    # stagetx <- stagetx %>% dplyr::rename(dexseq_stageR=transcript, feature_id=txID)
    # stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    # srestx <- full_join(stagetx, srestx, by=c("feature_id"="feature_id"))
    # restx <- restx %>% dplyr::rename(`dexseq_stageR`=transcript, feature_id=txID)

    write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

    # write.table(sresgene, paste0(outdir, sprintf("/results/stager_salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    # write.table(srestx, paste0(outdir, sprintf("/results/stager_salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
} 
### run kallisto counts

rm(genename)
rm(salmoncnt)
rm(cnt)
rm(d)
rm(resgene)
rm(restx)
rm(sresgene)
rm(srestx)

resgene <- read_tsv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)))

# sresgene <- read_tsv(paste0(outdir, sprintf("/results/stager_kal_res_gene_%s_%s.txt", con1, con2)))
# srestx <- read_tsv(paste0(outdir, sprintf("/results/stager_kal_res_tx_%s_%s.txt", con1, con2)))


if (!any(grepl("dexseq", colnames(resgene)))) {
    print("Running DEXSeq kallisto counts")
    gtf <- "/nfs/data/references/ensembl98_GRCh38/Homo_sapiens.GRCh38.98.gtf"
    tx2gene <- import_gtf(gtf_file = gtf)
    tx2gene$transcript_id_ver <- paste0(tx2gene$transcript_id,".", tx2gene$transcript_version)
    tx2gene <- move_columns_to_front(df = tx2gene, 
                                    columns = c("transcript_id_ver", "gene_id"))

    rfiles <- Sys.glob(paste0(outdir, "/kallisto_out/*/abundance.tsv"))
    names(rfiles) <- gsub(".*/","",gsub("/abundance.tsv","",rfiles))
    txi <- tximport(rfiles, type="kallisto", txOut=TRUE, tx2gene = tx2gene, ignoreAfterBar = TRUE)
    rsemtpm <- txi$counts
    rgenename <- read.csv(paste0(outdir, "/results/kal_count.csv"), sep="\t")
    rsemtpm <- apply(rsemtpm, 2, function(x){x+1})

    rsemcnt <- data.frame(gene_id=rgenename$gene_id, feature_id=rgenename$feature_id, rsemtpm)
    rsemcnt <- rsemcnt[rowSums(rsemtpm)>10,]
    colnames(rsemcnt) <- lapply(colnames(rsemcnt), function(x){gsub(".fq", "", gsub(".isoforms.results", "", x))}) %>% unlist

    rsemcnt <- rsemcnt %>% dplyr::filter(gene_id!="?")
    rsemcnt <- rsemcnt %>% dplyr::select("gene_id", "feature_id", meta1$sample_id)

    d <- dmDSdata(counts=rsemcnt, samples=meta1)

    dxd <- DEXSeqDataSet(countData=round(as.matrix(DRIMSeq::counts(d)[,-c(1:2)])), 
                        sampleData=DRIMSeq::samples(d),
                        design=~sample+exon+group:exon,
                        featureID=counts(d)$feature_id,
                        groupID=counts(d)$gene_id)

    dxd <- estimateSizeFactors(dxd)
    dxd <- estimateDispersions(dxd, BPPARAM = biocpar, quiet=T)
    dxd <- testForDEU(dxd, BPPARAM = biocpar, reducedModel = ~sample+exon)
    dxd = estimateExonFoldChanges( dxd, fitExpToVar="group", BPPARAM = biocpar)
    dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
    qval <- perGeneQValue(dxr)
    dxr <- as.data.frame(dxr)
    qval[qval==0] <- min(qval[qval !=0]) * 0.1
    dxr.g <- data.frame(gene=names(qval), qval)
    # dxr$pvalue[dxr$pvalue==0] <- min(dxr$pvalue[dxr$pvalue!=0]) * 0.1


    columns <- c("featureID","groupID","pvalue", "padj") 
    dxr <- as.data.frame(dxr[,columns])

    pConfirmation <- matrix(dxr$pvalue,ncol=1)
    dimnames(pConfirmation) <- list(c(dxr$featureID),c("transcript"))
    pScreen <- qval
    tx2gene <- counts(d)[,c("feature_id", "gene_id")]


    stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene) 
    stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=T) 
    suppressWarnings({ dex.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE) }) 



    #resgene <- inner_join(resgene, dxr %>% dplyr::select(groupID, pvalue), by=c("feature_id"="groupID"))
    dxr.g <- dxr.g %>% unique()
    resgene <- full_join(resgene, dxr.g, by=c("feature_id"="gene"))
    resgene <- resgene %>% dplyr::rename(`dexseq`=qval)

    # stagegene <- dex.padj %>% dplyr::select(geneID, gene) %>% unique()
    # stagegene <- stagegene %>% dplyr::rename(dexseq_stageR=gene, feature_id=geneID)
    # sresgene <- full_join(stagegene, sresgene, by=c("feature_id"="feature_id"))
    # #resgene <- resgene %>% dplyr::rename(`dexseq_stageR`=gene, feature_id=geneID)


    dxr.tx <- dxr %>% dplyr::select(featureID, padj) %>% dplyr::rename(dexseq=padj, feature_id=featureID)
    dxr.tx$feature_id <- lapply(dxr.tx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, dxr.tx, by=c("feature_id"="feature_id"))

    # stagetx <- dex.padj %>% dplyr::select(txID, transcript) %>% unique()
    # stagetx <- stagetx %>% dplyr::rename(dexseq_stageR=transcript, feature_id=txID)
    # stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    # srestx <- full_join(stagetx, srestx, by=c("feature_id"="feature_id"))
    # restx <- restx %>% dplyr::rename(`dexseq_stageR`=transcript, feature_id=txID)

    write.table(resgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

    # write.table(sresgene, paste0(outdir, sprintf("/results/stager_kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    # write.table(srestx, paste0(outdir, sprintf("/results/stager_kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}

rm(genename)
rm(salmoncnt)
rm(cnt)
rm(d)
rm(resgene)
rm(restx)
rm(sresgene)
rm(srestx)

resgene <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)))

# sresgene <- read_tsv(paste0(outdir, sprintf("/results/stager_rsem_res_gene_%s_%s.txt", con1, con2)))
# srestx <- read_tsv(paste0(outdir, sprintf("/results/stager_rsem_res_tx_%s_%s.txt", con1, con2)))

########## RSEM #########
if (!any(grepl("dexseq", colnames(resgene)))) {
    files <- Sys.glob(paste0(outdir, "/rsem_out/*/*.isoforms.results"))
    names(files) <- gsub(".*/","",gsub("/*.isoforms.results","",files))
    genename  <- read.csv(files[1], sep="\t")
    txi <- tximport(files, type="rsem", txOut=TRUE)
    rsemtpm <- txi$counts

    rsemcnt <- data.frame(gene_id=genename$gene_id, feature_id=genename$transcript_id, rsemtpm)
    rsemcnt <- rsemcnt[rowSums(rsemtpm)>10,]
    colnames(rsemcnt) <- lapply(colnames(rsemcnt), function(x){gsub(".fq", "", gsub(".isoforms.results", "", x))}) %>% unlist

    rsemcnt <- rsemcnt %>% dplyr::filter(gene_id!="?")
    rsemcnt <- rsemcnt %>% dplyr::select("gene_id", "feature_id", meta1$sample_id)

    d <- dmDSdata(counts=rsemcnt, samples=meta1)

    dxd <- DEXSeqDataSet(countData=round(as.matrix(DRIMSeq::counts(d)[,-c(1:2)])), 
                        sampleData=DRIMSeq::samples(d),
                        design=~sample+exon+group:exon,
                        featureID=counts(d)$feature_id,
                        groupID=counts(d)$gene_id)

    dxd <- estimateSizeFactors(dxd)
    dxd <- estimateDispersions(dxd, BPPARAM = biocpar, quiet=T)
    dxd <- testForDEU(dxd, BPPARAM = biocpar, reducedModel = ~sample+exon)
    dxd = estimateExonFoldChanges( dxd, fitExpToVar="group", BPPARAM = biocpar)
    dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
    qval <- perGeneQValue(dxr)
    dxr <- as.data.frame(dxr)
    qval[qval==0] <- min(qval[qval !=0]) * 0.1
    dxr.g <- data.frame(gene=names(qval), qval)
    # dxr$pvalue[dxr$pvalue==0] <- min(dxr$pvalue[dxr$pvalue!=0]) * 0.1


    columns <- c("featureID","groupID","pvalue", "padj") 
    dxr <- as.data.frame(dxr[,columns])

    pConfirmation <- matrix(dxr$pvalue,ncol=1)
    dimnames(pConfirmation) <- list(c(dxr$featureID),c("transcript"))
    pScreen <- qval
    tx2gene <- counts(d)[,c("feature_id", "gene_id")]


    stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene) 
    stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=T) 
    suppressWarnings({ dex.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=FALSE) }) 

    dxr.g <- dxr.g %>% unique()
    resgene <- full_join(resgene, dxr.g, by=c("feature_id"="gene"))
    resgene <- resgene %>% dplyr::rename(`dexseq`=qval)

    # stagegene <- dex.padj %>% dplyr::select(geneID, gene) %>% unique()
    # stagegene <- stagegene %>% dplyr::rename(dexseq_stageR=gene, feature_id=geneID)
    # sresgene <- full_join(stagegene, sresgene, by=c("feature_id"="feature_id"))
    # #resgene <- resgene %>% dplyr::rename(`dexseq_stageR`=gene, feature_id=geneID)


    dxr.tx <- dxr %>% dplyr::select(featureID, padj) %>% dplyr::rename(dexseq=padj, feature_id=featureID)
    dxr.tx$feature_id <- lapply(dxr.tx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, dxr.tx, by=c("feature_id"="feature_id"))

    # stagetx <- dex.padj %>% dplyr::select(txID, transcript) %>% unique()
    # stagetx <- stagetx %>% dplyr::rename(dexseq_stageR=transcript, feature_id=txID)
    # stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    # srestx <- full_join(stagetx, srestx, by=c("feature_id"="feature_id"))
    
    # restx <- restx %>% dplyr::rename(`drimseq_stageR`=transcript, feature_id=txID)

    write.table(resgene, paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

    # write.table(sresgene, paste0(outdir, sprintf("/results/stager_rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    # write.table(srestx, paste0(outdir, sprintf("/results/stager_rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}
