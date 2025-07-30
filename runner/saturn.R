library(tidyverse)
library(satuRn)
library(SummarizedExperiment)
library(tximport)
library(edgeR)
suppressMessages(library(DTUrtle))

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
path <- args[2]
meta <- args[3]
con1 <- args[4]
con2 <- args[5]
con <- c(con1, con2)
print(con)

resgene <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)))

# sresgene <- read_tsv(paste0(outdir, sprintf("/results/stager_salmon_res_gene_%s_%s.txt", con1, con2)))
# srestx <- read_tsv(paste0(outdir, sprintf("/results/stager_salmon_res_tx_%s_%s.txt", con1, con2)))

meta1 = read.csv(meta, sep="\t")
row.names(meta1) <- meta1$sample_id

meta1 <- meta1 %>% dplyr::filter(group %in% con)
meta1$sample_id <- lapply(meta1$sample_id, function(x){gsub(".sra", "", x)}) %>% unlist

if (!any(grepl("saturn", colnames(resgene)))) {
    print("Run saturn on salmon counts")
    genename  <- read.csv(paste0(outdir, "/results/salmon_count.csv"))
    txinfo <- genename %>% dplyr::select(gene_id, feature_id) %>% dplyr::rename(isoform_id = "feature_id")
    row.names(txinfo) <- txinfo$isoform_id
    files <- Sys.glob(paste0(outdir, "/salmon_out/*/quant.sf"))
    names(files) <- gsub(".*/","",gsub("/quant.sf","",files))
    txi <- tximport(files, type="salmon", txOut=TRUE)
    salmoncnt <- txi$counts

    salmoncnt <- salmoncnt[,meta1$sample_id]

    print("filter counts")
    filter_edgeR <- filterByExpr(salmoncnt,
                        design = NULL,
                        group = meta1$group,
                        lib.size = NULL,
                        min.count = 1,
                        min.total.count = 10,
                        large.n = 4,
                        min.prop = 0.7
                    ) # more stringen

    salmoncnt <- salmoncnt[filter_edgeR,]
    txinfo <- txinfo[which(txinfo$isoform_id %in% row.names(salmoncnt)),]
    txinfo <- subset(txinfo, 
                 duplicated(txinfo$gene_id) | duplicated(txinfo$gene_id, fromLast = TRUE))
    salmoncnt <- salmoncnt[which(row.names(salmoncnt) %in% txinfo$isoform_id), ]
    txinfo <- txinfo[match(row.names(salmoncnt), txinfo$isoform_id),]

    print("make exp")
    meta1$sample_name <- meta1$sample_id
    row.names(meta1) <- NULL
    sumExp <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = salmoncnt),
        colData = meta1,
        rowData = txinfo
    )

    metadata(sumExp)$formula <- ~ 0 + as.factor(colData(sumExp)$group)

    print("fitting DTU")
    sumExp <- satuRn::fitDTU(
                object = sumExp,
                formula = ~ 0 + group,
                parallel = FALSE,
                BPPARAM = BiocParallel::bpparam(),
                verbose = TRUE
            )

    group <- as.factor(meta1$group)
    design <- model.matrix(~0+group)
    colnames(design) <- levels(group)
    L <- limma::makeContrasts(Contrast1=sprintf("%s-%s", con1, con2), levels =design)

    print("testing DTU")
    sumExp <- satuRn::testDTU(
                object = sumExp,
                contrasts = L,
                sort = FALSE
            )

    print("joining results")
    res <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]
    res$isoform_id <- row.names(res)

    res <- full_join(res, txinfo, by=c("isoform_id"))
    res$isoform_id <- lapply(res$isoform_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist()
    
    res.g <- res %>% dplyr::select(gene_id, regular_FDR) %>% group_by(gene_id) %>% summarise(saturn=min(regular_FDR))

    resgene <- full_join(resgene, res.g, by=c("feature_id"="gene_id"))
    
    res.tx <- res %>% dplyr::select(isoform_id, regular_FDR) %>% dplyr::rename(saturn= regular_FDR)
    res.tx$isoform_id <- lapply(res.tx$isoform_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, res.tx, by=c("feature_id"="isoform_id"))

    write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

    # # transcript level p-values from satuRn
    # pvals <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]$empirical_pval

    #     # compute gene level q-values
    # geneID <- factor(rowData(sumExp)$gene_id)
    # geneSplit <- split(seq(along = geneID), geneID)
    # pGene <- sapply(geneSplit, function(i) min(pvals[i]))
    # pGene[is.na(pGene)] <- 1
    # theta <- unique(sort(pGene))

    # # gene-level significance testing
    # q <- DEXSeq:::perGeneQValueExact(pGene, theta, geneSplit) 
    # qScreen <- rep(NA_real_, length(pGene))
    # qScreen <- q[match(pGene, theta)]
    # qScreen <- pmin(1, qScreen)
    # names(qScreen) <- names(geneSplit)

    # # prepare stageR input
    # tx2gene <- as.data.frame(rowData(sumExp)[c("isoform_id", "gene_id")])
    # colnames(tx2gene) <- c("transcript", "gene")

    # pConfirmation <- matrix(matrix(pvals),
    #     ncol = 1,
    #     dimnames = list(rownames(tx2gene), "transcript")
    # )

    # # create a stageRTx object
    # stageRObj <- stageR::stageRTx(
    #     pScreen = qScreen,
    #     pConfirmation = pConfirmation,
    #     pScreenAdjusted = TRUE,
    #     tx2gene = tx2gene
    # )

    # # perform the two-stage testing procedure
    # stageRObj <- stageR::stageWiseAdjustment(
    #     object = stageRObj,
    #     method = "dtu",
    #     alpha = 0.05,
    #     allowNA = TRUE
    # )

    # # retrieves the adjusted p-values from the stageRTx object
    # padj <- stageR::getAdjustedPValues(stageRObj,
    #     order = TRUE,
    #     onlySignificantGenes = FALSE
    # )

    # stagegene <- padj %>% dplyr::select(geneID, gene) %>% unique()
    # stagegene <- stagegene %>% dplyr::rename(saturn_stageR=gene, feature_id=geneID)
    # sresgene <- full_join(stagegene, sresgene, by=c("feature_id"))

    # stagetx <- padj %>% dplyr::select(txID, transcript) %>% unique()
    # stagetx <- stagetx %>% dplyr::rename(saturn_stageR=transcript, feature_id=txID)
    # stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    # srestx <- full_join(stagetx, srestx, by=c("feature_id"))

    # write.table(sresgene, paste0(outdir, sprintf("/results/stager_salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    # # write.table(srestx, paste0(outdir, sprintf("/results/stager_salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}

rm(resgene)
rm(restx)
rm(skresgene)
rm(skrestx)

resgene <- read_tsv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)))

# skresgene <- read_tsv(paste0(outdir, sprintf("/results/stager_kal_res_gene_%s_%s.txt", con1, con2)))
# skrestx <- read_tsv(paste0(outdir, sprintf("/results/stager_kal_res_tx_%s_%s.txt", con1, con2)))

checkcount <- read_tsv(paste0(outdir, sprintf("/results/kal_count.csv")))
checkcount <- checkcount %>% dplyr::select(-feature_id, -gene_id)
colSums(checkcount)

if (!any(grepl("saturn", colnames(resgene)))) {
    print("Run saturn on kallisto counts")

    genename  <- read.csv(paste0(outdir, "/results/kal_count.csv"), sep="\t")
    txinfo <- genename %>% dplyr::select(gene_id, feature_id) %>% dplyr::rename(isoform_id = "feature_id")
    row.names(txinfo) <- txinfo$isoform_id
    print(meta1$sample_id) 

    gtf <- "/nfs/data/references/ensembl98_GRCh38/Homo_sapiens.GRCh38.98.gtf"
    tx2gene <- import_gtf(gtf_file = gtf)
    tx2gene$transcript_id_ver <- paste0(tx2gene$transcript_id,".", tx2gene$transcript_version)
    tx2gene <- move_columns_to_front(df = tx2gene, 
                                    columns = c("transcript_id_ver", "gene_id"))

    files <- paste0(outdir, "/kallisto_out/", meta1$sample_id,"/abundance.tsv")
    names(files) <- gsub(".*/","",gsub("/abundance.tsv","",files))
    txi <- tximport(files, type="kallisto", txOut=TRUE, tx2gene = tx2gene, ignoreAfterBar = TRUE)
    salmoncnt <- txi$counts

    salmoncnt <- salmoncnt[,meta1$sample_id]

    ##filter counts
    filter_edgeR <- filterByExpr(salmoncnt,
                        design = NULL,
                        group = meta1$group,
                        lib.size = NULL,
                        min.count = 1,
                        min.total.count = 10,
                        large.n = 4,
                        min.prop = 0.7
                    ) # more stringen

    salmoncnt <- salmoncnt[filter_edgeR,]
    txinfo <- txinfo[which(txinfo$isoform_id %in% row.names(salmoncnt)),]
    txinfo <- subset(txinfo, 
                 duplicated(txinfo$gene_id) | duplicated(txinfo$gene_id, fromLast = TRUE))
    salmoncnt <- salmoncnt[which(row.names(salmoncnt) %in% txinfo$isoform_id), ]
    txinfo <- txinfo[match(row.names(salmoncnt), txinfo$isoform_id),]

    meta1$sample_name <- meta1$sample_id
    row.names(meta1) <- NULL
    sumExp <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = salmoncnt),
        colData = meta1,
        rowData = txinfo
    )

    metadata(sumExp)$formula <- ~ 0 + as.factor(colData(sumExp)$group)

    sumExp <- satuRn::fitDTU(
                object = sumExp,
                formula = ~ 0 + group,
                parallel = FALSE,
                BPPARAM = BiocParallel::bpparam(),
                verbose = TRUE
            )

    group <- as.factor(meta1$group)
    design <- model.matrix(~0+group)
    colnames(design) <- levels(group)
    L <- limma::makeContrasts(Contrast1=sprintf("%s-%s", con1, con2), levels =design)

    sumExp <- satuRn::testDTU(
                object = sumExp,
                contrasts = L,
                sort = FALSE
            )

    res <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]
    res$isoform_id <- row.names(res)

    res <- full_join(res, txinfo, by=c("isoform_id"))
    res$isoform_id <- lapply(res$isoform_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist()
    
    res.g <- res %>% dplyr::select(gene_id, regular_FDR) %>% group_by(gene_id) %>% summarise(saturn=min(regular_FDR)) %>% ungroup 

    print("kallisto saturn")
    res.g %>% head
    
    resgene <- full_join(resgene, res.g, by=c("feature_id"="gene_id"))
    
    res.tx <- res %>% dplyr::select(isoform_id, regular_FDR) %>% rename(regular_FDR="saturn")
    res.tx$isoform_id <- lapply(res.tx$isoform_id, function(x){strsplit(x,"[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, res.tx, by=c("feature_id"="isoform_id"))

    write.table(resgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

    # transcript level p-values from satuRn
    # pvals <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]$empirical_pval

    #     # compute gene level q-values
    # geneID <- factor(rowData(sumExp)$gene_id)
    # geneSplit <- split(seq(along = geneID), geneID)
    # pGene <- sapply(geneSplit, function(i) min(pvals[i]))
    # pGene[is.na(pGene)] <- 1
    # theta <- unique(sort(pGene))

    # # gene-level significance testing
    # q <- DEXSeq:::perGeneQValueExact(pGene, theta, geneSplit) 
    # qScreen <- rep(NA_real_, length(pGene))
    # qScreen <- q[match(pGene, theta)]
    # qScreen <- pmin(1, qScreen)
    # names(qScreen) <- names(geneSplit)

    # # prepare stageR input
    # tx2gene <- as.data.frame(rowData(sumExp)[c("isoform_id", "gene_id")])
    # colnames(tx2gene) <- c("transcript", "gene")

    # pConfirmation <- matrix(matrix(pvals),
    #     ncol = 1,
    #     dimnames = list(rownames(tx2gene), "transcript")
    # )

    # # create a stageRTx object
    # stageRObj <- stageR::stageRTx(
    #     pScreen = qScreen,
    #     pConfirmation = pConfirmation,
    #     pScreenAdjusted = TRUE,
    #     tx2gene = tx2gene
    # )

    # # perform the two-stage testing procedure
    # stageRObj <- stageR::stageWiseAdjustment(
    #     object = stageRObj,
    #     method = "dtu",
    #     alpha = 0.05,
    #     allowNA = TRUE
    # )

    # # retrieves the adjusted p-values from the stageRTx object
    # padj <- stageR::getAdjustedPValues(stageRObj,
    #     order = TRUE,
    #     onlySignificantGenes = FALSE
    # )

    # stagegene <- padj %>% dplyr::select(geneID, gene) %>% unique()
    # stagegene <- stagegene %>% dplyr::rename(saturn_stageR=gene, feature_id=geneID)
    # skresgene <- full_join(stagegene, skresgene, by=c("feature_id"))

    # stagetx <- padj %>% dplyr::select(txID, transcript) %>% unique()
    # stagetx <- stagetx %>% dplyr::rename(saturn_stageR=transcript, feature_id=txID)
    # stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    # skrestx <- full_join(stagetx, skrestx, by=c("feature_id"))

    # write.table(skresgene, paste0(outdir, sprintf("/results/stager_kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    # write.table(skrestx, paste0(outdir, sprintf("/results/stager_kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

}

rm(resgene)
rm(restx)
rm(skresgene)
rm(skrestx)

resgene <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)))

# skresgene <- read_tsv(paste0(outdir, sprintf("/results/stager_rsem_res_gene_%s_%s.txt", con1, con2)))
# skrestx <- read_tsv(paste0(outdir, sprintf("/results/stager_rsem_res_tx_%s_%s.txt", con1, con2)))

if (!any(grepl("saturn", colnames(resgene)))) {
    print("Run saturn on rsem counts")
    
    
    files <- Sys.glob(paste0(outdir, "/rsem_out/*/*.isoforms.results"))
    names(files) <- gsub(".*/","",gsub("/*.isoforms.results","",files))
    genename  <- read.csv(files[1], sep="\t")
    txinfo <- genename %>% dplyr::select(gene_id, transcript_id) %>% dplyr::rename(isoform_id = "transcript_id")
    row.names(txinfo) <- txinfo$isoform_id
    txi <- tximport(files, type="rsem", txOut=TRUE)
    salmoncnt <- txi$counts

    salmoncnt <- salmoncnt[,meta1$sample_id]

    ##filter counts
    filter_edgeR <- filterByExpr(salmoncnt,
                        design = NULL,
                        group = meta1$group,
                        lib.size = NULL,
                        min.count = 1,
                        min.total.count = 10,
                        large.n = 4,
                        min.prop = 0.7
                    ) # more stringen

    salmoncnt <- salmoncnt[filter_edgeR,]
    txinfo <- txinfo[which(txinfo$isoform_id %in% row.names(salmoncnt)),]
    txinfo <- subset(txinfo, 
                 duplicated(txinfo$gene_id) | duplicated(txinfo$gene_id, fromLast = TRUE))
    salmoncnt <- salmoncnt[which(row.names(salmoncnt) %in% txinfo$isoform_id), ]
    txinfo <- txinfo[match(row.names(salmoncnt), txinfo$isoform_id),]

    meta1$sample_name <- meta1$sample_id
    row.names(meta1) <- NULL
    sumExp <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = salmoncnt),
        colData = meta1,
        rowData = txinfo
    )

    metadata(sumExp)$formula <- ~ 0 + as.factor(colData(sumExp)$group)

    sumExp <- satuRn::fitDTU(
                object = sumExp,
                formula = ~ 0 + group,
                parallel = FALSE,
                BPPARAM = BiocParallel::bpparam(),
                verbose = TRUE
            )

    group <- as.factor(meta1$group)
    design <- model.matrix(~0+group)
    colnames(design) <- levels(group)
    L <- limma::makeContrasts(Contrast1=sprintf("%s-%s", con1, con2), levels =design)

    sumExp <- satuRn::testDTU(
                object = sumExp,
                contrasts = L,
                sort = FALSE
            )

    res <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]
    res$isoform_id <- row.names(res)
    res <- full_join(res, txinfo, by=c("isoform_id"))
    res$isoform_id <- lapply(res$isoform_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist()
    
    res.g <- res %>% dplyr::select(gene_id, regular_FDR) %>% group_by(gene_id) %>% summarise(saturn=min(regular_FDR)) 
    resgene <- full_join(resgene, res.g, by=c("feature_id"="gene_id"))
    
    res.tx <- res %>% dplyr::select(isoform_id, regular_FDR) %>% rename(regular_FDR="saturn")
    res.tx$isoform_id <- lapply(res.tx$isoform_id, function(x){strsplit(x,"[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, res.tx, by=c("feature_id"="isoform_id"))

    write.table(resgene, paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

    # transcript level p-values from satuRn
    pvals <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]$empirical_pval

        # compute gene level q-values
    # geneID <- factor(rowData(sumExp)$gene_id)
    # geneSplit <- split(seq(along = geneID), geneID)
    # pGene <- sapply(geneSplit, function(i) min(pvals[i]))
    # pGene[is.na(pGene)] <- 1
    # theta <- unique(sort(pGene))

    # # gene-level significance testing
    # q <- DEXSeq:::perGeneQValueExact(pGene, theta, geneSplit) 
    # qScreen <- rep(NA_real_, length(pGene))
    # qScreen <- q[match(pGene, theta)]
    # qScreen <- pmin(1, qScreen)
    # names(qScreen) <- names(geneSplit)

    # # prepare stageR input
    # tx2gene <- as.data.frame(rowData(sumExp)[c("isoform_id", "gene_id")])
    # colnames(tx2gene) <- c("transcript", "gene")

    # pConfirmation <- matrix(matrix(pvals),
    #     ncol = 1,
    #     dimnames = list(rownames(tx2gene), "transcript")
    # )

    # # create a stageRTx object
    # stageRObj <- stageR::stageRTx(
    #     pScreen = qScreen,
    #     pConfirmation = pConfirmation,
    #     pScreenAdjusted = TRUE,
    #     tx2gene = tx2gene
    # )

    # # perform the two-stage testing procedure
    # stageRObj <- stageR::stageWiseAdjustment(
    #     object = stageRObj,
    #     method = "dtu",
    #     alpha = 0.05,
    #     allowNA = TRUE
    # )

    # # retrieves the adjusted p-values from the stageRTx object
    # padj <- stageR::getAdjustedPValues(stageRObj,
    #     order = TRUE,
    #     onlySignificantGenes = FALSE
    # )

    # stagegene <- padj %>% dplyr::select(geneID, gene) %>% unique()
    # stagegene <- stagegene %>% dplyr::rename(saturn_stageR=gene, feature_id=geneID)
    # skresgene <- full_join(stagegene, skresgene, by=c("feature_id"))

    # stagetx <- padj %>% dplyr::select(txID, transcript) %>% unique()
    # stagetx <- stagetx %>% dplyr::rename(saturn_stageR=transcript, feature_id=txID)
    # stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    # skrestx <- full_join(stagetx, skrestx, by=c("feature_id"))

    # write.table(skresgene, paste0(outdir, sprintf("/results/stager_rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    # write.table(skrestx, paste0(outdir, sprintf("/results/stager_rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

}


