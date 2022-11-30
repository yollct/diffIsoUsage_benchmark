library(tidyverse)
library(satuRn)
library(SummarizedExperiment)
library(tximport)
library(edgeR)

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

sresgene <- read_tsv(paste0(outdir, sprintf("/results/stager_salmon_res_gene_%s_%s.txt", con1, con2)))
srestx <- read_tsv(paste0(outdir, sprintf("/results/stager_salmon_res_tx_%s_%s.txt", con1, con2)))

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
                        min.total.count = 1,
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
                diagplot1 = TRUE,
                sort = FALSE
            )

    print("joining results")
    res <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]
    res$isoform_id <- row.names(res)

    res <- inner_join(res, txinfo, by=c("isoform_id"))
    res$isoform_id <- lapply(res$isoform_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist()
    
    res.g <- res %>% dplyr::select(gene_id, regular_FDR) %>% group_by(gene_id) %>% summarise(saturn=min(regular_FDR))

    resgene <- inner_join(resgene, res.g, by=c("feature_id"="gene_id"))
    
    res.tx <- res %>% dplyr::select(isoform_id, regular_FDR) %>% dplyr::rename(saturn= regular_FDR)
    res.tx$isoform_id <- lapply(res.tx$isoform_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- inner_join(restx, res.tx, by=c("feature_id"="isoform_id"))

    write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

    # transcript level p-values from satuRn
    pvals <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]$empirical_pval

        # compute gene level q-values
    geneID <- factor(rowData(sumExp)$gene_id)
    geneSplit <- split(seq(along = geneID), geneID)
    pGene <- sapply(geneSplit, function(i) min(pvals[i]))
    pGene[is.na(pGene)] <- 1
    theta <- unique(sort(pGene))

    # gene-level significance testing
    q <- DEXSeq:::perGeneQValueExact(pGene, theta, geneSplit) 
    qScreen <- rep(NA_real_, length(pGene))
    qScreen <- q[match(pGene, theta)]
    qScreen <- pmin(1, qScreen)
    names(qScreen) <- names(geneSplit)

    # prepare stageR input
    tx2gene <- as.data.frame(rowData(sumExp)[c("isoform_id", "gene_id")])
    colnames(tx2gene) <- c("transcript", "gene")

    pConfirmation <- matrix(matrix(pvals),
        ncol = 1,
        dimnames = list(rownames(tx2gene), "transcript")
    )

    # create a stageRTx object
    stageRObj <- stageR::stageRTx(
        pScreen = qScreen,
        pConfirmation = pConfirmation,
        pScreenAdjusted = TRUE,
        tx2gene = tx2gene
    )

    # perform the two-stage testing procedure
    stageRObj <- stageR::stageWiseAdjustment(
        object = stageRObj,
        method = "dtu",
        alpha = 0.05,
        allowNA = TRUE
    )

    # retrieves the adjusted p-values from the stageRTx object
    padj <- stageR::getAdjustedPValues(stageRObj,
        order = TRUE,
        onlySignificantGenes = FALSE
    )

    stagegene <- padj %>% dplyr::select(geneID, gene) %>% unique()
    stagegene <- stagegene %>% dplyr::rename(saturn_stageR=gene, feature_id=geneID)
    sresgene <- inner_join(stagegene, sresgene, by=c("feature_id"))

    stagetx <- padj %>% dplyr::select(txID, transcript) %>% unique()
    stagetx <- stagetx %>% dplyr::rename(saturn_stageR=transcript, feature_id=txID)
    stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    srestx <- inner_join(stagetx, srestx, by=c("feature_id"))

    write.table(sresggene, paste0(outdir, sprintf("/results/stager_salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(srestx, paste0(outdir, sprintf("/results/stager_salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}

rm(resgene)
rm(restx)
rm(sresgene)
rm(srestx)

resgene <- read_tsv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)))

skresgene <- read_tsv(paste0(outdir, sprintf("/results/stager_kal_res_gene_%s_%s.txt", con1, con2)))
skrestx <- read_tsv(paste0(outdir, sprintf("/results/stager_kal_res_tx_%s_%s.txt", con1, con2)))

if (!any(grepl("saturn", colnames(resgene)))) {
    print("Run saturn on kallisto counts")
    genename  <- read.csv(paste0(outdir, "/results/kal_count.csv"), sep="\t")
    txinfo <- genename %>% dplyr::select(gene_id, feature_id) %>% dplyr::rename(isoform_id = "feature_id")
    row.names(txinfo) <- txinfo$isoform_id
    files <- Sys.glob(paste0(outdir, "/kallisto_out/*/abundance.tsv"))
    names(files) <- gsub(".*/","",gsub("/abundance.tsv","",files))
    txi <- tximport(files, type="kallisto", txOut=TRUE)
    salmoncnt <- txi$counts

    salmoncnt <- salmoncnt[,meta1$sample_id]

    ##filter counts
    filter_edgeR <- filterByExpr(salmoncnt,
                        design = NULL,
                        group = meta1$group,
                        lib.size = NULL,
                        min.count = 1,
                        min.total.count = 1,
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
                diagplot1 = TRUE,
                sort = FALSE
            )

    res <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]
    res$isoform_id <- row.names(res)

    res <- inner_join(res, txinfo, by=c("isoform_id"))
    res$isoform_id <- lapply(res$isoform_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist()
    
    res.g <- res %>% dplyr::select(gene_id, regular_FDR) %>% group_by(gene_id) %>% summarise(saturn=min(regular_FDR)) %>% ungroup 
    resgene <- inner_join(resgene, res.g, by=c("feature_id"="gene_id"))
    
    res.tx <- res %>% dplyr::select(isoform_id, regular_FDR) %>% rename(regular_FDR="saturn")
    res.tx$isoform_id <- lapply(res.tx$isoform_id, function(x){strsplit(x,"[.]")[[1]][1]}) %>% unlist
    restx <- inner_join(restx, res.tx, by=c("feature_id"="isoform_id"))

    write.table(resgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

    # transcript level p-values from satuRn
    pvals <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]$empirical_pval

        # compute gene level q-values
    geneID <- factor(rowData(sumExp)$gene_id)
    geneSplit <- split(seq(along = geneID), geneID)
    pGene <- sapply(geneSplit, function(i) min(pvals[i]))
    pGene[is.na(pGene)] <- 1
    theta <- unique(sort(pGene))

    # gene-level significance testing
    q <- DEXSeq:::perGeneQValueExact(pGene, theta, geneSplit) 
    qScreen <- rep(NA_real_, length(pGene))
    qScreen <- q[match(pGene, theta)]
    qScreen <- pmin(1, qScreen)
    names(qScreen) <- names(geneSplit)

    # prepare stageR input
    tx2gene <- as.data.frame(rowData(sumExp)[c("isoform_id", "gene_id")])
    colnames(tx2gene) <- c("transcript", "gene")

    pConfirmation <- matrix(matrix(pvals),
        ncol = 1,
        dimnames = list(rownames(tx2gene), "transcript")
    )

    # create a stageRTx object
    stageRObj <- stageR::stageRTx(
        pScreen = qScreen,
        pConfirmation = pConfirmation,
        pScreenAdjusted = TRUE,
        tx2gene = tx2gene
    )

    # perform the two-stage testing procedure
    stageRObj <- stageR::stageWiseAdjustment(
        object = stageRObj,
        method = "dtu",
        alpha = 0.05,
        allowNA = TRUE
    )

    # retrieves the adjusted p-values from the stageRTx object
    padj <- stageR::getAdjustedPValues(stageRObj,
        order = TRUE,
        onlySignificantGenes = FALSE
    )

    stagegene <- padj %>% dplyr::select(geneID, gene) %>% unique()
    stagegene <- stagegene %>% dplyr::rename(saturn_stageR=gene, feature_id=geneID)
    skresgene <- inner_join(stagegene, skresgene, by=c("feature_id"))

    stagetx <- padj %>% dplyr::select(txID, transcript) %>% unique()
    stagetx <- stagetx %>% dplyr::rename(saturn_stageR=transcript, feature_id=txID)
    stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    skrestx <- inner_join(stagetx, skrestx, by=c("feature_id"))

    write.table(skresgene, paste0(outdir, sprintf("/results/stager_kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(skrestx, paste0(outdir, sprintf("/results/stager_kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

}

rm(resgene)
rm(restx)
rm(skresgene)
rm(skrestx)

resgene <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)))

skresgene <- read_tsv(paste0(outdir, sprintf("/results/stager_rsem_res_gene_%s_%s.txt", con1, con2)))
skrestx <- read_tsv(paste0(outdir, sprintf("/results/stager_rsem_res_tx_%s_%s.txt", con1, con2)))

if (!any(grepl("saturn", colnames(resgene)))) {
    print("Run saturn on rsem counts")
    
    txinfo <- genename %>% dplyr::select(gene_id, feature_id) %>% dplyr::rename(isoform_id = "feature_id")
    row.names(txinfo) <- txinfo$isoform_id
    files <- Sys.glob(paste0(outputdir, "/rsem_out/*/*.isoforms.resuts"))
    names(files) <- gsub(".*/","",gsub("/*.isoforms.results","",files))
    genename  <- read.csv(files[1], sep="\t")
    txi <- tximport(files, type="rsem", txOut=TRUE)
    salmoncnt <- txi$counts

    salmoncnt <- salmoncnt[,meta1$sample_id]

    ##filter counts
    filter_edgeR <- filterByExpr(salmoncnt,
                        design = NULL,
                        group = meta1$group,
                        lib.size = NULL,
                        min.count = 1,
                        min.total.count = 1,
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
                diagplot1 = TRUE,
                sort = FALSE
            )

    res <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]
    res$isoform_id <- row.names(res)
    res <- inner_join(res, txinfo, by=c("isoform_id"))
    res$isoform_id <- lapply(res$isoform_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist()
    
    res.g <- res %>% dplyr::select(gene_id, regular_FDR) %>% group_by(gene_id) %>% summarise(saturn=min(regular_FDR)) 
    resgene <- inner_join(resgene, res.g, by=c("feature_id"="gene_id"))
    
    res.tx <- res %>% dplyr::select(isoform_id, regular_FDR) %>% rename(regular_FDR="saturn")
    res.tx$isoform_id <- lapply(res.tx$isoform_id, function(x){strsplit(x,"[.]")[[1]][1]}) %>% unlist
    restx <- inner_join(restx, res.tx, by=c("feature_id"="isoform_id"))

    write.table(resgene, paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

    # transcript level p-values from satuRn
    pvals <- rowData(sumExp)[[sprintf("fitDTUResult_%s-%s", con1, con2)]]$empirical_pval

        # compute gene level q-values
    geneID <- factor(rowData(sumExp)$gene_id)
    geneSplit <- split(seq(along = geneID), geneID)
    pGene <- sapply(geneSplit, function(i) min(pvals[i]))
    pGene[is.na(pGene)] <- 1
    theta <- unique(sort(pGene))

    # gene-level significance testing
    q <- DEXSeq:::perGeneQValueExact(pGene, theta, geneSplit) 
    qScreen <- rep(NA_real_, length(pGene))
    qScreen <- q[match(pGene, theta)]
    qScreen <- pmin(1, qScreen)
    names(qScreen) <- names(geneSplit)

    # prepare stageR input
    tx2gene <- as.data.frame(rowData(sumExp)[c("isoform_id", "gene_id")])
    colnames(tx2gene) <- c("transcript", "gene")

    pConfirmation <- matrix(matrix(pvals),
        ncol = 1,
        dimnames = list(rownames(tx2gene), "transcript")
    )

    # create a stageRTx object
    stageRObj <- stageR::stageRTx(
        pScreen = qScreen,
        pConfirmation = pConfirmation,
        pScreenAdjusted = TRUE,
        tx2gene = tx2gene
    )

    # perform the two-stage testing procedure
    stageRObj <- stageR::stageWiseAdjustment(
        object = stageRObj,
        method = "dtu",
        alpha = 0.05,
        allowNA = TRUE
    )

    # retrieves the adjusted p-values from the stageRTx object
    padj <- stageR::getAdjustedPValues(stageRObj,
        order = TRUE,
        onlySignificantGenes = FALSE
    )

    stagegene <- padj %>% dplyr::select(geneID, gene) %>% unique()
    stagegene <- stagegene %>% dplyr::rename(saturn_stageR=gene, feature_id=geneID)
    skresgene <- inner_join(stagegene, skresgene, by=c("feature_id"))

    stagetx <- padj %>% dplyr::select(txID, transcript) %>% unique()
    stagetx <- stagetx %>% dplyr::rename(saturn_stageR=transcript, feature_id=txID)
    stagetx$feature_id <- lapply(stagetx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    skrestx <- inner_join(stagetx, skrestx, by=c("feature_id"))

    write.table(skresggene, paste0(outdir, sprintf("/results/stager_rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(skrestx, paste0(outdir, sprintf("/results/stager_rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

}

