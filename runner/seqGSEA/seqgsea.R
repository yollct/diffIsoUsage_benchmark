#install.packages("/nfs/home/students/chit/SeqGSEA_1.36.0.tar.gz", repos=NULL, type="source")
library(SeqGSEA)
library(DESeq2)
library(tidyverse)
library(doParallel)

cl <- makeCluster(20) # specify 2 cores to be used in computing
registerDoParallel(cl) # parallel backend registration
clusterEvalQ(cl, .libPaths("/nfs/home/students/chit/R/x86_64-pc-linux-gnu-library/4.2"))

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
path <- args[2]
meta <- args[3]
index <- args[4]
con1 <- args[5]
con2 <- args[6]
metadata <- read_tsv(meta)

ctrl.pattern <- con1# file name starting with "SC"
case.pattern <- con2 # file name starting with "SN"

 # gene set file and type
geneset.file <- paste0(index,"/this_gtf.exon")
data.dir <- paste0(outdir, "/results/htseq_exon")
output.prefix <- paste0(outdir, "/results/seqgsea")
geneID.type <- "ensembl"

case.files <- dir(data.dir, pattern=case.pattern, full.names = TRUE)
control.files <- dir(data.dir, pattern=ctrl.pattern, full.names = TRUE)
print("hier")
print(case.pattern)
print(control.files)

# analysis parameters
nCores <- 2
perm.times <- 100 # >= 1000 recommended
DEonly <- FALSE
DEweight <- c(0.2, 0.5, 0.8) # a vector for different weights
integrationMethod <- "linear"

### check files

RCS <- loadExonCountData(case.files, control.files)
RCS <- exonTestability(RCS, cutoff=5)
geneTestable <- geneTestability(RCS)
RCS <- subsetByGenes(RCS, unique(geneID(RCS))[ geneTestable ])
# get gene IDs, which will be used in initialization of gene set
geneIDs <- unique(geneID(RCS))
# calculate DS NB statistics

RCS <- estiExonNBstat(RCS)
RCS <- estiGeneNBstat(RCS)
# calculate DS NB statistics on the permutation data sets
permuteMat <- genpermuteMat(RCS, times=perm.times)
RCS <- DSpermute4GSEA(RCS, permuteMat)

# get gene read counts
geneCounts <- getGeneCount(RCS)
# calculate DE NB statistics
label <- label(RCS)
dds <-runDESeq(geneCounts, label)
DEGres <- DENBStat4GSEA(dds)

library(parallel)
cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, library(nlme))
clusterEvalQ(cl, library(DESeq2))
clusterExport(cl, ls())

times <- ncol(permuteMat)
n_gene <- nrow(counts(dds))
nbstat <- parLapply(cl=cl, X=1:times, fun=function(i) {
    newlabel <- as.factor(permuteMat[,i])
    dds <- DESeqDataSetFromMatrix(geneCounts, DataFrame(newlabel), ~ newlabel)
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    colA <- colData(dds)$newlabel == levels(colData(dds)$newlabel)[2]
    colB <- colData(dds)$newlabel == levels(colData(dds)$newlabel)[1]
    countsA <- counts(dds)[, colA]
    countsB <- counts(dds)[, colB]
    sizeFactorsA <- sizeFactors(dds)[colA]
    sizeFactorsB <- sizeFactors(dds)[colB]

    dispsA <- pmax(dispersions(dds), 1e-8)
    dispsB <- pmax(dispersions(dds), 1e-8)
    musA <- rowMeans(t(t(countsA)/sizeFactorsA))
    musB <- rowMeans(t(t(countsB)/sizeFactorsB))
    VarsA <- musA * sum(1/sizeFactorsA) / sum(colA)^2 + dispsA * musA^2 / sum(colA)
    VarsB <- musB * sum(1/sizeFactorsB) / sum(colB)^2 + dispsB * musB^2 / sum(colB)
    DEGresPerm <- data.frame(id = rownames(counts(dds)),
                    baseMeanA = musA, VarA = VarsA,
                    baseMeanB = musB, VarB = VarsB,
                    NBstat = (musA - musB) ^ 2 / (VarsA + VarsB),
                    stringsAsFactors = FALSE)
   DEGresPerm$NBstat
})

print("permuteNBstatGene")
DEpermNBstat <- do.call("cbind", nbstat)

DEscore.normFac <- normFactor(DEpermNBstat)
DEscore <- scoreNormalization(DEGres$NBstat, DEscore.normFac)
DEscore.perm <- scoreNormalization(DEpermNBstat, DEscore.normFac)

print("DS score normalization")
DSscore.normFac <- normFactor(RCS@permute_NBstat_gene)
DSscore <- scoreNormalization(RCS@featureData_gene$NBstat, DSscore.normFac)
DSscore.perm <- scoreNormalization(RCS@permute_NBstat_gene, DSscore.normFac)

RCS <- DSpermutePval(RCS, permuteMat)
rcs_seqgene  <- DSresultGeneTable(RCS)
RCS <- DSpermutePval(RCS, permuteMat)
seqgene <- rcs_seqgene %>% dplyr::select(geneID, padjust) %>% dplyr::rename("feature_id"="geneID", "seqGSEA"="padjust")

#plotGeneScore(gene.score, gene.score.perm, pdf=paste(output.prefix,".DSScore.pdf",sep=""), main="Splicing")

resgene <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)))
kalresgene <- read_tsv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)))
rsresgene <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)))
#resgene <- inner_join(resgene, dxr %>% dplyr::select(groupID, pvalue), by=c("feature_id"="groupID"))

# seqgene <- separate_rows(seqgene, feature_id, seqGSEA)
# seqgene %>% head

print("analysed seqGSEA")
print(dim(seqgene))
print(seqgene %>% head)

if (!any(grepl("seqGSEA", colnames(resgene)))) { 
    resgene <- full_join(resgene, seqgene, by="feature_id")
} 
if (!any(grepl("seqGSEA", colnames(kalresgene)))) { 
    kalresgene <- full_join(kalresgene, seqgene, by="feature_id")
}
if (!any(grepl("seqGSEA", colnames(rsresgene)))) { 
    rsresgene <- full_join(rsresgene, seqgene, by="feature_id")
}
#resgene <- resgene %>% dplyr::rename(`dexseq_stageR`=gene, feature_id=geneID)


write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
write.table(kalresgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
write.table(rsresgene, paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
