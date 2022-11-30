suppressMessages(library(IsoformSwitchAnalyzeR))
suppressMessages(library(tidyverse))
suppressMessages(library(AnnotationDbi))

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
path <- args[2]
meta <- args[3]
gtf <- args[4]

meta1 = read.csv(paste0(path,"/simulated_reads/sim_rep_info.txt"), sep="\t")
design <- data.frame(sampleID=meta1$rep_id,
                    condition=meta1$group)

### read salmon files
quant <- importIsoformExpression(parentDir = paste0(path, "/results/salmon_out/"))

##make tximeta
gtf1 <- rtracklayer::import(gtf)
gtf_df=as.data.frame(gtf1)


checkid = gtf_df$transcript_id[1]

for (x in 2:nrow(gtf_df)){
if (is.na(checkid)) {
  checkid <- gtf_df$transcript_id[x]
} else {
  break
}
}

# if (!grepl(".",checkid)){
#   print("make version")
gtf_df <- gtf_df %>% mutate(transcript_ver = paste0(transcript_id,".",transcript_version))

# } else {
#   gtf_df_filter <- gtf_df[gtf_df$transcript_id %in% quant$counts$isoform_id,]
#   rtracklayer::export(gtf_df_filter, '/nfs/proj/spycone/play_sample/filter_gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz')
# }

gtf_df_filter <- gtf_df[gtf_df$transcript_ver %in% quant$counts$isoform_id,]
quant$abundance <- quant$abundance[quant$abundance$isoform_id %in% gtf_df_filter$transcript_ver, ]
quant$counts <- quant$counts[quant$counts$isoform_id %in% gtf_df_filter$transcript_ver, ]
quant$length <- quant$length[quant$length$isoform_id %in% gtf_df_filter$transcript_ver, ]
#######
# fasta_set <- Biostrings::readDNAStringSet(filepath = "/nfs/proj/spycone/play_sample/gentrome_covid.fa.gz", format = 'fasta')

# names(fasta_set) <- unlist(lapply(names(fasta_set),function (x) {return(strsplit(x, " ")[[1]][1])}))
# quant$abundance <- quant$abundance[quant$abundance$isoform_id %in% names(fasta_set),]
# quant$counts <- quant$counts[quant$counts$isoform_id %in% names(fasta_set),]
# quant$length <- quant$length[quant$length$isoform_id %in% names(fasta_set),]


rtracklayer::export(gtf_df_filter %>% dplyr::select(-transcript_ver), paste0(path,"/annotations/filtered.gtf"), format="gtf")

aSwitchList <- importRdata(
  isoformCountMatrix   = quant$counts,
  isoformRepExpression = quant$abundance,
  designMatrix         = design,
  isoformExonAnnoation = paste0(path,"/annotations/filtered.gtf"),
  isoformNtFasta       = fasta,
  ignoreAfterPeriod = TRUE,
  ignoreAfterSpace = TRUE,
  ignoreAfterBar = TRUE,
  showProgress = TRUE
)

exampleSwitchListAnalyzed.g <- isoformSwitchTestDEXSeq(
      switchAnalyzeRlist = aSwitchList,
      reduceToSwitchingGenes=TRUE,
      dIFcutoff = 0,
      alpha=1
)


exampleSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
      switchAnalyzeRlist = aSwitchList,
      reduceToSwitchingGenes=FALSE,
      dIFcutoff = 0,
      alpha=1
)
######to do implement gene level 

genemapper <- gtf_df_filter %>% dplyr::select(gene_id, gene_name) %>% unique()
txmapper <- gtf_df_filter %>% dplyr::select(transcript_id, transcript_ver)
isa_dex_g <- exampleSwitchListAnalyzed.g$isoformFeatures %>% dplyr::select(gene_name, gene_switch_q_value)
isa_dex_g <- inner_join(isa_dex_g, genemapper %>% dplyr::select(gene_name, gene_id), by="gene_name")
isa_dex_g$gene_switch_q_value[isa_dex_g$gene_switch_q_value==0] <- min(isa_dex_g$gene_switch_q_value[isa_dex_g$gene_switch_q_value!=0]) * 0.1

isa_dex_tx <- exampleSwitchListAnalyzed$isoformFeatures %>% filter(isoform_switch_q_value < 1 & dIF>=0)
isa_dex_tx <- inner_join(isa_dex_tx, genemapper %>% dplyr::select(gene_name, gene_id), by="gene_name")
isa_dex_tx$gene_switch_q_value[isa_dex_tx$gene_switch_q_value==0] <- min(isa_dex_tx$gene_switch_q_value[isa_dex_tx$gene_switch_q_value!=0]) * 0.1

##### add the result to result table
resgene <- read_tsv(paste0(path, "/results/res_gene.txt"))
restx <- read_tsv(paste0(path,"/results/res_tx.csv"))

isa_dex_g <- isa_dex_g %>% dplyr::rename(isa_dexseq=gene_switch_q_value, feature_id=gene_id) %>% unique()
resgene <- inner_join(resgene, isa_dex_g, by="feature_id")

#resgene <- inner_join(resgene, dxr %>% dplyr::select(groupID, pvalue), by=c("feature_id"="groupID"))
isa_dex_tx <- inner_join(isa_dex_tx, txmapper, by=c("isoform_id"="transcript_id"))
isa_dex_tx <- isa_dex_tx %>% dplyr::select(transcript_ver, isoform_switch_q_value) %>% dplyr::rename(isa_dexseq=isoform_switch_q_value, feature_id=transcript_ver) %>% unique()
restx <- full_join(restx, isa_dex_tx, by=c("feature_id"="feature_id"))


write.table(resgene, paste0(path, "/results/res_gene.txt"), row.names = FALSE, sep="\t")
write.table(restx, paste0(path, "/results/res_tx.csv"), row.names=FALSE, sep="\t")

# pr_dexseq <- data.frame()
# tx_pr_dexseq <- data.frame()
# tx_gg_pr_dexseq <- data.frame()
# gg_pr_dexseq <- data.frame()
# geneobj_dexseq <- list()

# ei0.1g <- filter_eventim(type="gene")
# ei0.1t <- filter_eventim(type="isoform")
# i<-1
# for (co in c(0.01,0.05,0.1)){
#   for (dif in seq(0,1,0.01)){
#     for (evn in c(0, 1)){
#       try(isa_dex_tx <- exampleSwitchListAnalyzed$isoformFeatures %>% filter(isoform_switch_q_value < co & dIF>=dif))
#       isa_dex_tx_geneid <- AnnotationDbi::select(org.Hs.eg.db, keys=isa_dex_tx$gene_name, keytype="SYMBOL", columns="ENSEMBL")
#       if (evn==0){
#         pr <- overall_precision_recall(isa_dex_tx_geneid$ENSEMBL, groundtruth, IS_tx)
#         gg_pr <- gene_group_precision_recall(isa_dex_tx_geneid$ENSEMBL, groundtruth, IS_tx)
#         tx_pr <- overall_tx_precision_recall(isa_dex_tx$isoform_id, tx_groundtruth, tx_gg)
#         tx_gg_pr <- gene_group_tx_precision_recall(isa_dex_tx$isoform_id, tx_groundtruth, tx_gg)
#       } else if (evn==1){
#         pr <- overall_precision_recall(isa_dex_tx_geneid$ENSEMBL[isa_dex_tx_geneid$ENSEMBL %in% ei0.1g$gene], groundtruth, IS_tx)
#         gg_pr <- gene_group_precision_recall(isa_dex_tx_geneid$ENSEMBL[isa_dex_tx_geneid$ENSEMBL %in% ei0.1g$gene], groundtruth, IS_tx)
#         tx_pr <- overall_tx_precision_recall(isa_dex_tx$isoform_id[isa_dex_tx$isoform_id %in% ei0.1t$isoform], tx_groundtruth, tx_gg)
#         tx_gg_pr <- gene_group_tx_precision_recall(isa_dex_tx$isoform_id[isa_dex_tx$isoform_id %in% ei0.1t$isoform], tx_groundtruth, tx_gg)
#       }
#       #save(isa_dex_tx, file=paste0(data, "/results/isa_dex_tx.rda"))
      
#       if (co==0.05 & dif == 0){
#         geneobj_dexseq[[i]] <- isa_dex_tx_geneid$ENSEMBL
#         names(geneobj_dexseq[[i]]) <- sprintf("ISA (DEXSeq)_%s", dif)
#         i <- i+1
#       }
      
      
#       gg_pr$tool <- sprintf("ISA (DEXSeq)_%s", dif)
#       gg_pr$cutoff <- co
#       gg_pr$eventim <- ifelse(evn==1, "yes", "no")
#       gg_pr_dexseq <- rbind(gg_pr_dexseq, gg_pr)
#       pr$tool <- sprintf("ISA (DEXSeq)_%s", dif)
#       pr$cutoff <- co
#       pr$eventim <- ifelse(evn==1, "yes", "no")
#       pr_dexseq <- rbind(pr_dexseq, pr)
#       tx_pr$tool <- sprintf("ISA (DEXSeq)_%s", dif)
#       tx_pr$cutoff <- co
#       tx_pr$eventim <- ifelse(evn==1, "yes", "no")
#       tx_pr_dexseq <- rbind(tx_pr_dexseq, tx_pr)
#       tx_gg_pr$tool <- sprintf("ISA (DEXSeq)_%s", dif)
#       tx_gg_pr$cutoff <- co
#       tx_gg_pr$eventim <- ifelse(evn==1, "yes", "no")
#       tx_gg_pr_dexseq <- rbind(tx_gg_pr_dexseq, tx_gg_pr)
#       rm(isa_dex_tx)
#     }
#   }
# }

# exampleSwitchListAnalyzed <- isoformSwitchTestDRIMSeq(
#   switchAnalyzeRlist = aSwitchList,
#   testIntegration='isoform_only',
#   reduceToSwitchingGenes=TRUE,
#   dIFcutoff = 0,
#   alpha=1
# )

# pr_drimseq <- data.frame()
# gg_pr_drimseq <- data.frame()
# tx_pr_drimseq <- data.frame()
# tx_gg_pr_drimseq <- data.frame()
# geneobj_drimseq <- list()
# i<-1
# for (co in c(0.01,0.05,0.1)){
#   for (dif in seq(0,1,0.01)){
#     for (evn in c(0,1)){
#     try(isa_drim_tx <- exampleSwitchListAnalyzed$isoformFeatures %>% filter(isoform_switch_q_value < co & dIF>=dif))
#     isa_drim_tx_geneid <- AnnotationDbi::select(org.Hs.eg.db, keys=isa_drim_tx$gene_name, keytype="SYMBOL", columns="ENSEMBL")
#     if (evn==0){
#     #save(isa_drim_tx, file=paste0(data, "/results/isa_drim_tx.rda"))
#       pr <- overall_precision_recall(isa_drim_tx_geneid$ENSEMBL, groundtruth, IS_tx)
#       gg_pr <- gene_group_precision_recall(isa_drim_tx_geneid$ENSEMBL, groundtruth, IS_tx)
#       tx_pr <- overall_tx_precision_recall(isa_drim_tx$isoform_id, tx_groundtruth, tx_gg)
#       tx_gg_pr <- gene_group_tx_precision_recall(isa_drim_tx$isoform_id, tx_groundtruth, tx_gg)
#     } else if (evn==1){
#       pr <- overall_precision_recall(isa_drim_tx_geneid$ENSEMBL[isa_drim_tx_geneid$ENSEMBL %in% ei0.1g$gene], groundtruth, IS_tx)
#       gg_pr <- gene_group_precision_recall(isa_drim_tx_geneid$ENSEMBL[isa_drim_tx_geneid$ENSEMBL %in% ei0.1g$gene], groundtruth, IS_tx)
#       tx_pr <- overall_tx_precision_recall(isa_drim_tx$isoform_id[isa_drim_tx$isoform_id %in% ei0.1t$isoform], tx_groundtruth, tx_gg)
#       tx_gg_pr <- gene_group_tx_precision_recall(isa_drim_tx$isoform_id[isa_drim_tx$isoform_id %in% ei0.1t$isoform], tx_groundtruth, tx_gg)
#     }
#     #save(isa_dex_tx, file=paste0(data, "/results/isa_dex_tx.rda"))
#     pr$tool <- sprintf("ISA (DRIMSeq)_%s", dif)
#     pr$cutoff <- co
#     pr$eventim <- ifelse(evn==1, "yes", "no")
#     pr_drimseq <- rbind(pr_drimseq, pr)

#     if (co==0.05 & dif == 0){
#       geneobj_drimseq[[i]] <- isa_drim_tx_geneid$ENSEMBL
#       names(geneobj_drimseq[[i]]) <- sprintf("ISA (DRIMSeq)_%s", dif)
#       i <- i+1
#     }

#     gg_pr$tool <- sprintf("ISA (DRIMSeq)_%s", dif)
#     gg_pr$cutoff <- co
#     gg_pr$eventim <- ifelse(evn==1, "yes", "no")
#     gg_pr_drimseq <- rbind(gg_pr_drimseq, gg_pr)
#     tx_gg_pr$tool <- sprintf("ISA (DRIMSeq)_%s", dif)
#     tx_gg_pr$cutoff <- co
#     tx_gg_pr$eventim <- ifelse(evn==1, "yes", "no")
#     tx_gg_pr_drimseq <- rbind(tx_gg_pr_drimseq, tx_gg_pr)
#     tx_pr$tool <- sprintf("ISA (DRIMSeq)_%s", dif)
#     tx_pr$cutoff <- co
#     tx_pr$eventim <- ifelse(evn==1, "yes", "no")
#     tx_pr_drimseq <- rbind(tx_pr_drimseq, tx_pr)
#     }
# }
# }


# isa_pr <- rbind(pr_dexseq, pr_drimseq)
# write.table(isa_pr, paste0(path,"/results/isa_pr_df.csv"), row.names=F)
# isa_gg_pr <- rbind(gg_pr_dexseq, gg_pr_drimseq)
# write.table(isa_gg_pr, paste0(path,"/results/isa_gg_pr_df.csv"), row.names=F)
# isa_tx_pr <- rbind(tx_pr_dexseq, tx_pr_drimseq)
# write.table(isa_tx_pr, paste0(path,"/results/isa_tx_pr_df.csv"), row.names=F)
# isa_tx_gg_pr <- rbind(tx_gg_pr_dexseq, tx_gg_pr_drimseq)
# write.table(isa_tx_gg_pr, paste0(path,"/results/isa_tx_gg_pr_df.csv"), row.names=F)

# saveRDS(geneobj_dexseq, paste0(path,"/results/isa_dexseq_0.rds"))
# saveRDS(geneobj_drimseq, paste0(path,"/results/isa_drimseq_0.rds"))


# g<-switchPlot(
#   exampleSwitchListAnalyzed,
#   gene='SELENOP',
#   localTheme = theme_bw(base_size = 13) # making text sightly larger for vignette
# )
# ggsave(g, file=paste0(path, "/results/switchplot.png"))

# tx <- extractTopSwitches(exampleSwitchListAnalyzed, filterForConsequences = FALSE, n=NA, dIFcutoff=0.1,alpha=0.1, extractGenes = F)
# overall_tx_precision_recall(tx$isoform_id, tx_groundtruth, tx_gg)

# isa_dex_tx_geneid$ENSEMBL[!isa_dex_tx_geneid$ENSEMBL %in% groundtruth$geneid]

