library(limma)
library(edgeR)
suppressMessages(library(tximport))
suppressMessages(library(tidyverse))
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

print("in limma")
resgene <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)))

meta1 = read.csv(meta, sep="\t")
row.names(meta1) <- meta1$sample_id
meta1$sample_id <- lapply(meta1$sample_id, function(x){gsub(".sra", "",x)}) %>% unlist

if (!any(grepl("LimmaDS", colnames(restx)))) { 
    genename  <- read.csv(paste0(outdir, "/results/salmon_count.csv"))
    files <- Sys.glob(paste0(outdir, "/salmon_out/*/quant.sf"))
    names(files) <- gsub(".*/","",gsub("/quant.sf","",files))
    txi <- tximport(files, type="salmon", txOut=TRUE)
    salmontpm <- txi$counts

    salmoncnt <- data.frame(txi$counts)
    rownames(salmoncnt) <- genename$feature_id

    dge <- DGEList(counts=salmoncnt,genes=genename$feature_id)
    dge$genes$GeneID <- genename$gene_id

    A <- rowSums(dge$counts)
    dge <- dge[A>10, keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)

    design <- model.matrix(~ group, data = meta1)
    v <- voom(dge, design, plot=F)
    fit <- lmFit(v, design)

    ex <- diffSplice(fit, geneid="GeneID")
    res <- topSplice(ex, n=nrow(dge$genes), coef=2, test="t")
 
    res1 <- res %>% dplyr::select(GeneID, FDR) %>% group_by(GeneID)  %>% summarise(FDR=min(FDR)) %>% na.omit() %>% unique()
    limmares <- res1 %>% dplyr::rename("feature_id"="GeneID", "LimmaDS"="FDR")

    resgene$feature_id <- as.character(resgene$feature_id)
    resgene <- full_join(resgene, limmares, by="feature_id")

    limmatx <- res %>% dplyr::select(genes, FDR) %>% dplyr::filter(FDR < 0.05) %>% dplyr::rename("feature_id"="genes", "LimmaDS"="FDR")
    limmatx$feature_id <- lapply(limmatx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, limmatx, by="feature_id")

    write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
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

if (!any(grepl("LimmaDS", colnames(restx)))) {
    print("Running Limma kallisto counts")
    gtf <- "/nfs/data/references/ensembl98_GRCh38/Homo_sapiens.GRCh38.98.gtf"
    tx2gene <- import_gtf(gtf_file = gtf)
    tx2gene$transcript_id_ver <- paste0(tx2gene$transcript_id,".", tx2gene$transcript_version)
    tx2gene <- move_columns_to_front(df = tx2gene, 
                                    columns = c("transcript_id_ver", "gene_id"))
                                    
    rfiles <- Sys.glob(paste0(outdir, "/kallisto_out/*/abundance.tsv"))
    names(rfiles) <- gsub(".*/","",gsub("/abundance.tsv","",rfiles))
    txi <- tximport(rfiles, type="kallisto", txOut=TRUE, tx2gene = tx2gene, ignoreAfterBar = TRUE)

    rgenename <- read.csv(paste0(outdir, "/results/kal_count.csv"), sep="\t")

    cnt <- data.frame(txi$counts)
    rownames(cnt) <- rgenename$feature_id

    dge <- DGEList(counts=cnt,genes=rgenename$feature_id)
    dge$genes$GeneID <- rgenename$gene_id

    A <- rowSums(dge$counts)
    dge <- dge[A>10,, keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)

    design <- model.matrix(~ group, data = meta1)
    v <- voom(dge, design, plot=F)
    fit <- lmFit(v, design)

    ex <- diffSplice(fit, geneid="GeneID")
    res <- topSplice(ex, n=nrow(dge$genes), coef=2, test="t")

    res1 <- res %>% dplyr::select(GeneID, FDR) %>% group_by(GeneID) %>% summarise(FDR=min(FDR)) %>% na.omit() %>% unique() 
    limmares <- res1 %>% dplyr::rename("feature_id"="GeneID", "LimmaDS"="FDR")

    resgene$feature_id <- as.character(resgene$feature_id)
    resgene <- full_join(resgene, limmares, by="feature_id")

    limmatx <- res %>% dplyr::select(genes, FDR) %>% dplyr::filter(FDR < 0.05) %>% dplyr::rename("feature_id"="genes", "LimmaDS"="FDR")
    limmatx$feature_id <- lapply(limmatx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, limmatx, by="feature_id")

    write.table(resgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
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

########## RSEM #########
if (!any(grepl("LimmaDS", colnames(restx)))) {
    files <- Sys.glob(paste0(outdir, "/rsem_out/*/*.isoforms.results"))
    names(files) <- gsub(".*/","",gsub("/*.isoforms.results","",files))
    genename  <- read.csv(files[1], sep="\t")
    txi <- tximport(files, type="rsem", txOut=TRUE)

    cnt <- data.frame(txi$counts)
    rownames(cnt) <- genename$transcript_id

    dge <- DGEList(counts=cnt, genes=genename$transcript_id)
    dge$genes$GeneID <- genename$gene_id

    A <- rowSums(dge$counts)
    dge <- dge[A>10,, keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)

    design <- model.matrix(~ group, data = meta1)
    v <- voom(dge, design, plot=F)
    fit <- lmFit(v, design)

    ex <- diffSplice(fit, geneid="GeneID")
    res <- topSplice(ex, n=nrow(dge$genes), coef=2, test="t")

    res1 <- res %>% dplyr::select(GeneID, FDR) %>% group_by(GeneID) %>% summarise(FDR=min(FDR)) %>% na.omit() %>% unique() 
    limmares <- res1 %>% dplyr::rename("feature_id"="GeneID", "LimmaDS"="FDR")

    resgene$feature_id <- as.character(resgene$feature_id)
    resgene <- full_join(resgene, limmares, by="feature_id")


    limmatx <- res %>% dplyr::select(genes, FDR) %>% dplyr::filter(FDR < 0.05) %>% dplyr::rename("feature_id"="genes", "LimmaDS"="FDR")
    limmatx$feature_id <- lapply(limmatx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, limmatx, by="feature_id")

    write.table(resgene, paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
}
