suppressMessages(library(tidyverse))
suppressMessages(library(DTUrtle))
biocpar <- BiocParallel::MulticoreParam(4)
#or standard serial computation (only 1 core)
biocpar <- BiocParallel::SerialParam()

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
path <- args[2]
meta <- args[3]
gtf <- args[4]
con1 <- args[5]
con2 <- args[6]
con <- c(con1, con2)
print(con)

##### add the result to result table
resgene <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)))

meta1 = read.csv(meta, sep="\t")
row.names(meta1) <- meta1$sample

if (!any(grepl("dturtle", colnames(resgene)))) { 
    meta1 <- meta1 %>% dplyr::filter(group %in% con)

    # source(paste0(path,"/runner/groundtruth.R"))
    # source(paste0(path,"/runner/event_importance.R"))
    #the BiocParallel framework is used to parallelize the computations.
    #Using 4 cores:
    
    #multiple other options available for computational clusters.


    # #import gtf Annotation to get transcript to gene mapping
    tx2gene <- import_gtf(gtf_file = gtf)
    tx2gene$transcript_id_ver <- paste0(tx2gene$transcript_id,".", tx2gene$transcript_version)
    tx2gene <- move_columns_to_front(df = tx2gene, 
                                        columns = c("transcript_id_ver", "gene_id"))

    files <- Sys.glob(paste0(outdir, "/salmon_out/*/quant.sf"))
    names(files) <- gsub(".*/","",gsub("/quant.sf","",files))
    cts <- import_counts(files = files, type = "salmon")

    pd <- data.frame("id"=colnames(cts), "group"=meta1$group, 
                    stringsAsFactors = FALSE)

    dturtle <- run_drimseq(counts = cts, tx2gene = tx2gene, pd=pd, id_col = "id",
                        cond_col = "group", filtering_strategy = "bulk", 
                        BPPARAM = biocpar, subset_feature = row.names(cts) %in% tx2gene$transcript_id_ver)

    # ei0.1g <- filter_eventim(type="gene")
    # ei0.1t <- filter_eventim(type="isoform")

    dturtle_1 <- posthoc_and_stager(dturtle = dturtle, ofdr = 1, posthoc = 0)
    saveRDS(dturtle, paste0(outdir, "/results/salmon_dturtle_res.rds"))
    dturtle <- readRDS(paste0(outdir, "/results/salmon_dturtle_res.rds"))
    names(dturtle_1)
    dturtle_1$FDR_table %>% dim
    dturtle_1$sig_tx[1:10]


    dturtle_g <- data.frame(feature_id=dturtle_1$FDR_table$geneID, dturtle=dturtle_1$FDR_table$gene) %>% unique
    resgene <- inner_join(resgene, dturtle_g, by="feature_id")

    #resgene <- inner_join(resgene, dxr %>% dplyr::select(groupID, pvalue), by=c("feature_id"="groupID"))
    dturtle_tx <- data.frame(feature_id=dturtle_1$FDR_table$txID, dturtle=dturtle_1$FDR_table$transcript)
    dturtle_tx$feature_id <- lapply(dturtle_tx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, dturtle_tx, by=c("feature_id"="feature_id"))

    write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}

####### RSEM counts ###############
rm(resgene)
rm(restx)
rm(sresgene)
rm(srestx)
rm(pd)
rm(cts)

resgene <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)))

if (!any(grepl("dturtle", colnames(resgene)))) { 
    print("dturtle to rsem")
    meta1 <- meta1 %>% dplyr::filter(group %in% con)

    # source(paste0(path,"/runner/groundtruth.R"))
    # source(paste0(path,"/runner/event_importance.R"))
    #the BiocParallel framework is used to parallelize the computations.
    #Using 4 cores:
    
    #multiple other options available for computational clusters.


    # #import gtf Annotation to get transcript to gene mapping
    tx2gene <- import_gtf(gtf_file = gtf)
    tx2gene$transcript_id_ver <- paste0(tx2gene$transcript_id,".", tx2gene$transcript_version)
    tx2gene <- move_columns_to_front(df = tx2gene, 
                                        columns = c("transcript_id", "gene_id"))

    files <- Sys.glob(paste0(outdir, "/rsem_out/*/*.isoforms.results"))
    names(files) <- gsub(".*/","",gsub("/*.isoforms.results","",files))
    cts <- import_counts(files = files, type = "rsem")

    pd <- data.frame("id"=colnames(cts), "group"=meta1$group, 
                    stringsAsFactors = FALSE)

    dturtle <- run_drimseq(counts = cts, tx2gene = tx2gene, pd=pd, id_col = "id",
                        cond_col = "group", filtering_strategy = "bulk", 
                        BPPARAM = biocpar, subset_feature = row.names(cts) %in% tx2gene$transcript_id)

    # ei0.1g <- filter_eventim(type="gene")
    # ei0.1t <- filter_eventim(type="isoform")

    dturtle_1 <- posthoc_and_stager(dturtle = dturtle, ofdr = 1, posthoc = 0)
    saveRDS(dturtle, paste0(outdir, "/results/rsem_dturtle_res.rds"))
    dturtle <- readRDS(paste0(outdir, "/results/rsem_dturtle_res.rds"))
    names(dturtle_1)
    dturtle_1$FDR_table %>% dim
    dturtle_1$sig_tx[1:10]


    dturtle_g <- data.frame(feature_id=dturtle_1$FDR_table$geneID, dturtle=dturtle_1$FDR_table$gene) %>% unique
    resgene <- inner_join(resgene, dturtle_g, by="feature_id")

    #resgene <- inner_join(resgene, dxr %>% dplyr::select(groupID, pvalue), by=c("feature_id"="groupID"))
    dturtle_tx <- data.frame(feature_id=dturtle_1$FDR_table$txID, dturtle=dturtle_1$FDR_table$transcript)
    dturtle_tx$feature_id <- lapply(dturtle_tx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, dturtle_tx, by=c("feature_id"="feature_id"))

    write.table(resgene, paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}

rm(resgene)
rm(restx)
rm(sresgene)
rm(srestx)
rm(pd)
rm(cts)

resgene <- read_tsv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)))


if (!any(grepl("dturtle", colnames(resgene)))) { 
    #import gtf Annotation to get transcript to gene mapping
    tx2gene <- import_gtf(gtf_file = gtf)
    tx2gene$transcript_id_ver <- paste0(tx2gene$transcript_id,".", tx2gene$transcript_version)
    tx2gene <- move_columns_to_front(df = tx2gene, 
                                        columns = c("transcript_id_ver", "gene_id"))

    rfiles <- Sys.glob(paste0(outdir, "/kallisto_out/*/abundance.h5"))
    names(rfiles) <- gsub(".*/","",gsub("/abundance.h5","",rfiles))
    rcts <- import_counts(files = rfiles, type = "kallisto")

    rpd <- data.frame("id"=colnames(rcts), "group"=meta1$group, 
                    stringsAsFactors = FALSE)

    dturtle <- run_drimseq(counts = rcts, tx2gene = tx2gene, pd=rpd, id_col = "id",
                        cond_col = "group", filtering_strategy = "bulk", 
                        BPPARAM = biocpar, subset_feature = row.names(rcts) %in% tx2gene$transcript_id_ver)

    # ei0.1g <- filter_eventim(type="gene")
    # ei0.1t <- filter_eventim(type="isoform")

    dturtle_1 <- posthoc_and_stager(dturtle = dturtle, ofdr = 1, posthoc = 0)
    saveRDS(dturtle, paste0(outdir, "/results/kal_dturtle_res.rds"))
    dturtle <- readRDS(paste0(outdir, "/results/kal_dturtle_res.rds"))


    ##### add the result to result table
    resgene <- read_tsv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)))
    restx <- read_tsv(paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)))

    dturtle_g <- data.frame(feature_id=dturtle_1$FDR_table$geneID, dturtle=dturtle_1$FDR_table$gene) %>% unique
    resgene <- inner_join(resgene, dturtle_g, by="feature_id")

    #resgene <- inner_join(resgene, dxr %>% dplyr::select(groupID, pvalue), by=c("feature_id"="groupID"))
    dturtle_tx <- data.frame(feature_id=dturtle_1$FDR_table$txID, dturtle=dturtle_1$FDR_table$transcript)
    dturtle_tx$feature_id <- lapply(dturtle_tx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
    restx <- full_join(restx, dturtle_tx, by=c("feature_id"="feature_id"))

    write.table(resgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
    write.table(restx, paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")
}


# pr_dturtle <- data.frame()
# gg_pr_dturtle <- data.frame()
# tx_pr_dturtle <- data.frame()
# tx_gg_pr_dturtle <- data.frame()
# geneobj_dturtle <- list()
# i<-1
# for (co in c(0.01,0.05,0.1)){
#     for (evn in c(0,1)){
#     dif <- 0
#     dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = co)

#     if (evn==0){
#     #save(isa_drim_tx, file=paste0(data, "/results/isa_drim_tx.rda"))
#       pr <- overall_precision_recall(dturtle$sig_gene, groundtruth, IS_tx)
#       gg_pr <- gene_group_precision_recall(dturtle$sig_gene, groundtruth, IS_tx)
#       tx_pr <- overall_tx_precision_recall(dturtle$sig_tx, tx_groundtruth, tx_gg)
#       tx_gg_pr <- gene_group_tx_precision_recall(dturtle$sig_tx, tx_groundtruth, tx_gg)
#     } else if (evn==1){
#       pr <- overall_precision_recall(dturtle$sig_gene[dturtle$sig_gene %in% ei0.1g$gene], groundtruth, IS_tx)
#       gg_pr <- gene_group_precision_recall(dturtle$sig_gene[dturtle$sig_gene %in% ei0.1g$gene], groundtruth, IS_tx)
#       tx_pr <- overall_tx_precision_recall(dturtle$sig_tx[dturtle$sig_tx %in% ei0.1t$isoform], tx_groundtruth, tx_gg)
#       tx_gg_pr <- gene_group_tx_precision_recall(dturtle$sig_tx[dturtle$sig_tx %in% ei0.1g$isoform], tx_groundtruth, tx_gg)
#     }
#     #save(isa_dex_tx, file=paste0(data, "/results/isa_dex_tx.rda"))
#     pr$tool <- sprintf("dturtle_%s", dif)
#     pr$cutoff <- co
#     pr$eventim <- ifelse(evn==1, "yes", "no")
#     pr_dturtle <- rbind(pr_dturtle, pr)

#     if (co==0.05 & dif == 0){
#       geneobj_dturtle[[i]] <- dturtle$sig_gene
#       names(geneobj_dturtle[[i]]) <- sprintf("dturtle_%s", dif)
#       i <- i+1
#     }

#     gg_pr$tool <- sprintf("dturtle_%s", dif)
#     gg_pr$cutoff <- co
#     gg_pr$eventim <- ifelse(evn==1, "yes", "no")
#     gg_pr_dturtle <- rbind(gg_pr_dturtle, gg_pr)
#     tx_gg_pr$tool <- sprintf("dturtle_%s", dif)
#     tx_gg_pr$cutoff <- co
#     tx_gg_pr$eventim <- ifelse(evn==1, "yes", "no")
#     tx_gg_pr_dturtle <- rbind(tx_gg_pr_dturtle, tx_gg_pr)
#     tx_pr$tool <- sprintf("dturtle_%s", dif)
#     tx_pr$cutoff <- co
#     tx_pr$eventim <- ifelse(evn==1, "yes", "no")
#     tx_pr_dturtle <- rbind(tx_pr_dturtle, tx_pr)
#     }
# }


# write.table(pr_dturtle, paste0(path,"/results/dtu_pr_df.csv"), row.names=F)
# write.table(gg_pr, paste0(path,"/results/dtu_gg_pr_df.csv"), row.names=F)
# write.table(tx_pr_dturtle, paste0(path,"/results/dtu_tx_pr_df.csv"), row.names=F)
# write.table(tx_gg_pr_dturtle, paste0(path,"/results/dtu_tx_gg_pr_df.csv"), row.names=F)

# saveRDS(geneobj_dturtle, paste0(path,"/results/dtu_0.rds"))