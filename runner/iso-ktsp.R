library(tidyverse)


args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
path <- args[2]
con1 <- args[3]
con2 <- args[4]
print(con1)
print(con2)

ktsp_res <- readLines(paste0(outdir, sprintf("/results/salmon_isoktsp_output_%s_%s.txt", con1, con2)))

resgene <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)))

if (!any(grepl("iso_ktsp", colnames(resgene)))) { 
    ktsp_gene <- c()
    gene_score <- c()
    ktsp_isoform <- c()
    for (k in ktsp_res){
        if (strsplit(k, "\t")[[1]][2]=="single_pair_performance"){
            ktsp_gene <- c(ktsp_gene, strsplit(strsplit(k, "\t")[[1]][3], ",")[[1]][1])
            ktsp_isoform <- c(ktsp_isoform, strsplit(strsplit(k, "\t")[[1]][3], ",")[[1]][2])
            ktsp_isoform <- c(ktsp_isoform, strsplit(strsplit(k, "\t")[[1]][4], ",")[[1]][2])
            gene_score <- c(gene_score, strsplit(strsplit(k, "\t")[[1]][9], "=")[[1]][2])
        }
    }
    iso_score <- rep(gene_score, each=2)

    gene_score <- as.numeric(gene_score)
    iso_score <- as.numeric(iso_score)
    ktsp_g <- data.frame(feature_id=ktsp_gene, `iso_ktsp`=sapply(gene_score,function(x){1-x}))
    ktsp_iso <- data.frame(feature_id=ktsp_isoform, `iso_ktsp`=sapply(iso_score, function(x){1-x}))

    ktsp_g %>% head

    resgene <- full_join(ktsp_g, resgene, by="feature_id")
    restx <- full_join(ktsp_iso, restx, by="feature_id")

    write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}

ktsp_res = readLines(paste0(outdir, sprintf("/results/kal_isoktsp_output_%s_%s.txt", con1, con2)))

rm(resgene)
rm(restx)
rm(ktsp_res)
ktsp_res <- readLines(paste0(outdir, sprintf("/results/kal_isoktsp_output_%s_%s.txt", con1, con2)))
resgene <- read_tsv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)))

if (!any(grepl("iso_ktsp", colnames(resgene)))) { 
    ktsp_gene <- c()
    gene_score <- c()
    ktsp_isoform <- c()
    for (k in ktsp_res){
        if (strsplit(k, "\t")[[1]][2]=="single_pair_performance"){
            ktsp_gene <- c(ktsp_gene, strsplit(strsplit(k, "\t")[[1]][3], ",")[[1]][1])
            ktsp_isoform <- c(ktsp_isoform, strsplit(strsplit(k, "\t")[[1]][3], ",")[[1]][2])
            ktsp_isoform <- c(ktsp_isoform, strsplit(strsplit(k, "\t")[[1]][4], ",")[[1]][2])
            gene_score <- c(gene_score, strsplit(strsplit(k, "\t")[[1]][9], "=")[[1]][2])
        }
    }
    iso_score <- rep(gene_score, each=2)

    gene_score <- as.numeric(gene_score)
    iso_score <- as.numeric(iso_score)
    ktsp_g <- data.frame(feature_id=ktsp_gene, `iso_ktsp`=sapply(gene_score,function(x){1-x}))
    ktsp_iso <- data.frame(feature_id=ktsp_isoform, `iso_ktsp`=sapply(iso_score, function(x){1-x}))




    resgene <- full_join(ktsp_g, resgene, by="feature_id")
    restx <- full_join(ktsp_iso, restx, by="feature_id")

    write.table(resgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}


rm(resgene)
rm(restx)
rm(ktsp_res)
ktsp_res <- readLines(paste0(outdir, sprintf("/results/rsem_isoktsp_output_%s_%s.txt", con1, con2)))
resgene <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)))

if (!any(grepl("iso_ktsp", colnames(resgene)))) { 
    ktsp_gene <- c()
    gene_score <- c()
    ktsp_isoform <- c()
    for (k in ktsp_res){
        if (strsplit(k, "\t")[[1]][2]=="single_pair_performance"){
            ktsp_gene <- c(ktsp_gene, strsplit(strsplit(k, "\t")[[1]][3], ",")[[1]][1])
            ktsp_isoform <- c(ktsp_isoform, strsplit(strsplit(k, "\t")[[1]][3], ",")[[1]][2])
            ktsp_isoform <- c(ktsp_isoform, strsplit(strsplit(k, "\t")[[1]][4], ",")[[1]][2])
            gene_score <- c(gene_score, strsplit(strsplit(k, "\t")[[1]][9], "=")[[1]][2])
        }
    }
    iso_score <- rep(gene_score, each=2)

    gene_score <- as.numeric(gene_score)

    iso_score <- as.numeric(iso_score)

    ktsp_g <- data.frame(feature_id=ktsp_gene, `iso_ktsp`=sapply(gene_score,function(x){1-x}))
    ktsp_iso <- data.frame(feature_id=ktsp_isoform, `iso_ktsp`=sapply(iso_score, function(x){1-x}))
    
    resgene <- full_join(ktsp_g, resgene, by="feature_id")
    restx <- full_join(ktsp_iso, restx, by="feature_id")

    write.table(resgene, paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}
# pr_dtu <- data.frame()
# gg_pr_dtu <- data.frame()
# tx_pr_dtu <- data.frame()
# tx_gg_pr_dtu <- data.frame()
# ei0.1g <- filter_eventim(type="gene")
# ei0.1t <- filter_eventim(type="isoform")

# for (ofdr in c(50, 100, 150)){
#     for (evn in c(0,1)){
#         if (evn==0){
#             pr <- overall_precision_recall(ktsp_gene[1:ofdr], groundtruth, IS_tx)
#             gg_pr <- gene_group_precision_recall(ktsp_gene[1:ofdr], groundtruth, IS_tx)
#             tx_pr <- overall_tx_precision_recall(ktsp_isoform[1:ofdr*2], tx_groundtruth, tx_gg)
#             tx_gg_pr <- gene_group_tx_precision_recall(ktsp_isoform[1:ofdr*2], tx_groundtruth, tx_gg)
#         } else if (evn==1){
#             pr <- overall_precision_recall(ktsp_gene[1:ofdr][ktsp_gene[1:ofdr] %in% ei0.1g$gene], groundtruth, IS_tx)
#             gg_pr <- gene_group_precision_recall(ktsp_gene[1:ofdr][ktsp_gene[1:ofdr] %in% ei0.1g$gene], groundtruth, IS_tx)
#             tx_pr <- overall_tx_precision_recall(ktsp_isoform[1:ofdr*2][ktsp_isoform[1:ofdr*2] %in% ei0.1t$isoform], tx_groundtruth, tx_gg)
#             tx_gg_pr <- gene_group_tx_precision_recall(ktsp_isoform[1:ofdr*2][ktsp_isoform[1:ofdr*2] %in% ei0.1t$isoform], tx_groundtruth, tx_gg)
#         }

#         pr$tool <- sprintf("iso-KTSP")
#         pr$cutoff <- ofdr
#         pr$eventim <- ifelse(evn==1, "yes", "no")
#         pr_dtu <- rbind(pr_dtu, pr)

#         gg_pr$tool <- sprintf("iso-KTSP")
#         gg_pr$cutoff <- ofdr
#         gg_pr$eventim <- ifelse(evn==1, "yes", "no")
#         gg_pr_dtu <- rbind(gg_pr_dtu, gg_pr)
#         tx_pr$tool <- sprintf("iso-KTSP")
#         tx_pr$cutoff <- ofdr
#         tx_pr$eventim <- ifelse(evn==1, "yes", "no")
#         tx_gg_pr$tool <- sprintf("iso-KTSP")
#         tx_gg_pr$cutoff <- ofdr
#         tx_gg_pr$eventim <- ifelse(evn==1, "yes", "no")
#         tx_pr_dtu <- rbind(tx_pr_dtu, tx_pr)
#         tx_gg_pr_dtu <- rbind(tx_gg_pr_dtu, tx_gg_pr)
#     }
# }

# write.table(pr_dtu, paste0(path,"/results/ktsp_pr_df.csv"), row.names=F)
# write.table(gg_pr_dtu, paste0(path,"/results/ktsp_gg_pr_df.csv"), row.names=F)
# write.table(tx_pr_dtu, paste0(path,"/results/ktsp_tx_pr_df.csv"), row.names=F)
# write.table(tx_gg_pr_dtu, paste0(path,"/results/ktsp_tx_gg_pr_df.csv"), row.names=F)
##compare the two tools