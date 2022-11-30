library(dplyr)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly=TRUE)

outdir <- args[1]

overall_precision_recall <- function(genelist, groundtruth, allgenelist){
  tp <- sum(genelist %in% groundtruth$geneid)
  fp <- sum(!genelist %in% groundtruth$geneid)
  #tn <- sum(!allgenelist$geneid %in% c(genelist, groundtruth$geneid))
  fn <- sum(!groundtruth$geneid %in% genelist)

  precision  <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  return(data.frame(precision=precision, recall=recall))
}

gene_group_precision_recall <- function(genelist, groundtruth, allgenelist){
  thisgenelist <- allgenelist[allgenelist$geneid %in% genelist,]
  pr <- c()
  re <- c()
  grp <- c()
  for (gp in c("2-4", "5-9", ">9")){
    tmp_gt <- groundtruth %>% dplyr::filter(gene_group==gp)
    tmp_all <- allgenelist %>% dplyr::filter(gene_group==gp)
    tmp_gene <- thisgenelist %>% dplyr::filter(gene_group==gp)

    if (dim(tmp_gene)[1]==0){
      precision <- NA
      recall <- NA
    }
    tp <- sum(tmp_gene$geneid %in% tmp_gt$geneid)
    fp <- sum(!tmp_gene$geneid %in% tmp_gt$geneid)
    tn <- sum(!tmp_all$geneid %in% c(tmp_gene$geneid, tmp_gt$geneid))
    fn <- sum(!tmp_gt$geneid %in% tmp_gene$geneid)
    
    precision  <- tp/(tp+fp)
    recall <- tp/(tp+fn)
    
    pr <- c(pr, precision)
    re <- c(re, recall)
    grp <- c(grp, gp)
  }

  res_group <- data.frame(precision=pr, recall=re, gene_group=grp)
  return(res_group)
}

#### get groundtruth (in gene level)
##### only genes with 2 DE isoforms are considered as isoform switched genes
sim_tx <- read.csv(paste0(outdir, "/simulated_reads/sim_tx_info.txt"), sep="\t")
sim_tx$geneid <- lapply(sim_tx$transcriptid, function(x){strsplit(str_replace(strsplit(x, split=" ")[[1]][4], "gene:", ""),"[.]")[[1]][1]}) %>% unlist
sim_tx$transcriptid <- lapply(sim_tx$transcriptid, function(x){strsplit(x, " ")[[1]][1]}) %>% unlist
sim_tx$DE.status.bin <- sim_tx$DEstatus.V1 + sim_tx$DEstatus.c2
sim_tx$fc <- pmax(sim_tx$foldchange.V1, sim_tx$foldchange.c2)
IS_tx <- sim_tx %>% group_by(geneid) %>% mutate(IS=sum(DE.status.bin), n_tx=n(), fc=max(fc)) %>% ungroup
IS_tx$geneid <- lapply(IS_tx$geneid, function(x){strsplit(strsplit(x, " ")[[1]][1], "[.]")[[1]][1]}) %>% unlist
IS_tx$genesymb <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                         keys=IS_tx$geneid, 
                                         keytype="ENSEMBL",
                                        column="SYMBOL",
                                        multiVals = "first")
IS_tx <- IS_tx %>% dplyr::mutate(gene_group = ifelse(n_tx>9, ">9", ifelse(n_tx>5, "5-9", ifelse(n_tx>2, "2-4", "1"))))

#groundtruth <- IS_tx %>% filter(IS>1)

truthtable_gene <- sim_tx %>% group_by(geneid) %>% summarise(IS=sum(DE.status.bin), n_tx=n(), fc=max(fc)) %>% dplyr::mutate(status=ifelse(IS>=2, 1, 0), diffiso_group=ifelse(IS<=1, 0, ifelse(IS>10, ">10", ifelse(IS>=6, "6-9", ifelse(IS>=3, "3-5", 2)))), fc_group=ifelse(fc==2, "2", ifelse(fc<5, "3-5", ifelse(fc<8, "6-8", ">8")))) 
truthtable_gene <- truthtable_gene %>% dplyr::rename(feature_id=geneid)
write.table(truthtable_gene, paste0(outdir, sprintf("/results/truthtable_gene.csv")), sep="\t", row.names = F)

truthtable_tx <- IS_tx %>% dplyr::mutate(status=ifelse(DE.status.bin==1,ifelse(IS>=2, 1, 0),0), diffiso_group=ifelse(IS<=1, 0, ifelse(IS>10, ">10", ifelse(IS>=6, "6-9", ifelse(IS>=3, "3-5", 2)))), fc_group=ifelse(fc==2, "2", ifelse(fc<5, "3-5", ifelse(fc<8, "6-8", ">8")))) 
truthtable_tx <- truthtable_tx %>% dplyr::rename(feature_id=transcriptid)
write.table(truthtable_tx, paste0(outdir, sprintf("/results/truthtable_tx.csv")), sep="\t", row.names = F)


# tx_gg <- inner_join(sim_tx, IS_tx %>% dplyr::select(geneid, gene_group), by="geneid")

# tx_groundtruth <- sim_tx %>% dplyr::filter(geneid %in% groundtruth$geneid & DE.status.bin>0) %>% inner_join(tx_gg, by="transcriptid")

# overall_tx_precision_recall <- function(txlist, groundtruth, allgenelist){
#   tp <- sum(txlist %in% groundtruth$transcriptid)
#   fp <- sum(!txlist %in% groundtruth$transcriptid)
#   #tn <- sum(!allgenelist$geneid %in% c(genelist, groundtruth$geneid))
#   fn <- sum(!groundtruth$transcriptid %in% txlist)

#   precision  <- tp/(tp+fp)
#   recall <- tp/(tp+fn)
#   return(data.frame(precision=precision, recall=recall))
# }

# gene_group_tx_precision_recall <- function(genelist, groundtruth, allgenelist){
#   thisgenelist <- allgenelist[allgenelist$transcriptid %in% genelist,]
#   pr <- c()
#   re <- c()
#   grp <- c()
#   for (gp in c("2-4", "5-9", ">9")){
#     tmp_gt <- groundtruth %>% dplyr::filter(gene_group==gp)
#     tmp_all <- allgenelist %>% dplyr::filter(gene_group==gp)
#     tmp_gene <- thisgenelist %>% dplyr::filter(gene_group==gp)

#     if (dim(tmp_gene)[1]==0){
#       precision <- NA
#       recall <- NA
#     }
#     tp <- sum(tmp_gene$transcriptid %in% tmp_gt$transcriptid)
#     fp <- sum(!tmp_gene$transcriptid %in% tmp_gt$transcriptid)
#     tn <- sum(!tmp_all$transcriptid %in% c(tmp_gene$transcriptid, tmp_gt$transcriptid))
#     fn <- sum(!tmp_gt$transcriptid %in% tmp_gene$transcriptid)
    
#     precision  <- tp/(tp+fp)
#     recall <- tp/(tp+fn)
    
#     pr <- c(pr, precision)
#     re <- c(re, recall)
#     grp <- c(grp, gp)
#   }

#   res_group <- data.frame(precision=pr, recall=re, gene_group=grp)
#   return(res_group)
# }

