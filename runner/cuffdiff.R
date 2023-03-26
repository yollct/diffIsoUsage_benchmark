library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
path <- args[2]
meta <- args[3]
con1 <- args[4]
con2 <- args[5]
gtfpath <- args[6]
#gtfpath <- "/nfs/data/references/ensembl98_GRCh38/Homo_sapiens.GRCh38.98.gtf"
con <- c(con1, con2)
print(con)

genediff <- read.csv(paste0(outdir, "/results/cuffdiff_results/gene_exp.diff"), sep="\t")
isodiff <- read.csv(paste0(outdir, "/results/cuffdiff_results/isoform_exp.diff"), sep="\t")
anno <- read.csv(paste0(outdir, "/results/cuffdiff_results/isoforms.fpkm_tracking"), sep="\t")
isodiff <- isodiff %>% full_join(anno, by=c("test_id"="tracking_id"))

resgene <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)))
kresgene <- read_tsv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)))
krestx <- read_tsv(paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)))
rsresgene <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)))
rsrestx <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)))

gtf <- rtracklayer::import(gtfpath)
gtf <- as.data.frame(gtf)
gtf <- gtf %>% dplyr::filter(type=="gene")
# sym2ens <- unique(data.frame(symb=gtf$gene_name, ensembl=gtf$gene_id))
# row.names(sym2ens) <- sym2ens$symb
#isodiff_id <- AnnotationDbi::select(org.Hs.eg.db, keys=isodiff$gene, keytype="SYMBOL", columns="ENSEMBL")
cuffdiff_tx <- data.frame(feature_id=anno$nearest_ref_id, cuffdiff=isodiff$q_value) %>% unique()

cuffdiff_g <- data.frame(feature_id=genediff$gene_id, cuffdiff=genediff$q_value) 

if (!any(grepl("cuffdiff", colnames(resgene)))) { 
    resgene <- full_join(resgene, cuffdiff_g, by="feature_id")
    restx <- full_join(restx, cuffdiff_tx, by="feature_id")
}
if (!any(grepl("cuffdiff", colnames(kresgene)))) { 
    kresgene <- full_join(kresgene, cuffdiff_g, by="feature_id")
    krestx <- full_join(krestx, cuffdiff_tx, by="feature_id")
}
if (!any(grepl("cuffdiff", colnames(rsresgene)))) {
    rsresgene <- full_join(rsresgene, cuffdiff_g, by="feature_id")
    rsrestx <- full_join(rsrestx, cuffdiff_tx, by="feature_id")
}

write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
write.table(restx, paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
write.table(kresgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
write.table(krestx, paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
write.table(rsresgene, paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
write.table(rsrestx, paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
# cuffidff_g <- data.frame()
# pr_cd <- data.frame()
# pr_gg_cd <- data.frame()
# tx_pr_cd <- data.frame()
# tx_gg_pr_cd <-data.frame()
# for (co in c(0.01,0.05, 0.1)){

#     sigdiff <- isodiff %>% dplyr::filter(q_value<co)
#     sigdiff_id <- AnnotationDbi::select(org.Hs.eg.db, keys=sigdiff$gene, keytype="SYMBOL", columns="ENSEMBL")
#     pr <- overall_precision_recall(sigdiff_id$ENSEMBL, groundtruth, IS_tx)
#     gg_pr <- gene_group_precision_recall(sigdiff_id$ENSEMBL, groundtruth, IS_tx)
#     tx_pr <- overall_tx_precision_recall(sigdiff$test_id, tx_groundtruth, tx_gg)
#     tx_gg_pr <- gene_group_tx_precision_recall(sigdiff$test_id, tx_groundtruth, tx_gg)
    
#     pr$cutoff <- co
#     pr_cd <- rbind(pr_cd, pr)
#     gg_pr$cutoff <- co
#     pr_gg_cd <- rbind(pr_gg_cd, gg_pr)
#     tx_pr$cutoff <- co
#     tx_pr_cd <- rbind(tx_pr_cd, tx_pr)
#     tx_gg_pr$cutoff <- co
#     tx_gg_pr_cd <- rbind(tx_gg_pr_cd, tx_gg_pr)
# }

# pr_cd$tool <- "Cuffdiff"
# pr_gg_cd$tool <- "Cuffdiff"
# tx_pr_cd$tool <- "Cuffdiff"
# tx_gg_pr_cd$tool <- "Cuffdiff"

# write.table(pr_cd, paste0(res,"/results/cd_pr_df.csv"), row.names=F)
# write.table(pr_gg_cd, paste0(res,"/results/cd_gg_pr_df.csv"), row.names=F)
# write.table(tx_pr_cd, paste0(res,"/results/cd_tx_pr_df.csv"), row.names=F)
# write.table(tx_gg_pr_cd, paste0(res,"/results/cd_tx_gg_pr_df.csv"), row.names=F)
