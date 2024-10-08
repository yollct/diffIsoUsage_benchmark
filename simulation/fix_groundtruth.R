library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
path <- "/nfs/home/students/chit/is_benchmark"
outdir <- sprintf("/nfs/scratch/chit/new_simulations/%s/", dir)

sim_g <- read.csv(paste0(outdir, "results/geneinfo.csv"), sep="\t")
sim_g$tmp <- paste0(as.character(sim_g$DTE),as.character(sim_g$IS),as.character(sim_g$DTU))
sim_g <- sim_g %>% mutate(events=ifelse(tmp=="100", "DTE", ifelse(tmp=="010", "IS", ifelse(tmp=="001", "DTU", "none")))) %>% dplyr::select(events, feature_id)
sim_g %>% head
sim_tx <- read.csv(paste0(outdir, "results/isoinfo.csv"), sep="\t")
sim_tx %>% head
sims <- left_join(sim_tx, sim_g, by=c("gene_id"="feature_id"))



# sim_tx$geneid <- lapply(sim_tx$transcriptid, function(x){strsplit(str_replace(strsplit(x, split=" ")[[1]][4], "gene:", ""),"[.]")[[1]][1]}) %>% unlist
#sim_tx$transcript_id <- lapply(sim_tx$transcriptid, function(x){strsplit(x, " ")[[1]][1]}) %>% unlist
#sim_tx$DE.status.bin <- sim_tx$DEstatus.V1 + sim_tx$DEstatus.c2
#sim_tx$fc <- pmax(sim_tx$foldchange.V1, sim_tx$foldchange.c2)
sim_tx1 <- sims %>% group_by(gene_id) %>% mutate(IS=sum(DE.status.bin), n_tx=n(), fc=max(fc)) %>% ungroup


sim_tx1 <- sim_tx1 %>% dplyr::mutate(gene_group = ifelse(n_tx>9, ">9", ifelse(n_tx>5, "5-8", ifelse(n_tx>2, "2-4", "1"))))


#groundtruth <- IS_tx %>% filter(IS>1)

truthtable_gene <- sim_tx1 %>% group_by(gene_id) %>% summarise(events=max(events), IS=sum(DE.status.bin), n_tx=n(), fc=max(fc)) %>% dplyr::mutate(status=ifelse(IS>=1, 1, 0), diffiso_group=ifelse(IS<=1, 0, ifelse(IS>10, ">10", ifelse(IS>=6, "6-9", ifelse(IS>=3, "3-5", 2)))), gene_group = ifelse(n_tx>9, ">9", ifelse(n_tx>=5, "5-9", ifelse(n_tx>=2, "2-4", "1"))))  %>% dplyr::rename(fc_group="fc")
truthtable_gene <- truthtable_gene %>% dplyr::rename(feature_id=gene_id)
truthtable_gene %>% head
write.table(truthtable_gene, paste0(outdir, sprintf("/results/truthtable_gene.csv")), sep="\t", row.names = F)

truthtable_tx <- sim_tx %>% group_by(gene_id) %>% dplyr::mutate(IS=sum(DE.status.bin)) %>% ungroup() %>% dplyr::mutate(status=ifelse(DE.status.bin==1,ifelse(IS>=2, 1, 0),0), diffiso_group=ifelse(IS<=1, 0, ifelse(IS>10, ">10", ifelse(IS>=6, "6-9", ifelse(IS>=3, "3-5", 2))))) 
truthtable_tx <- truthtable_tx %>% dplyr::rename(feature_id=transcript_id)
write.table(truthtable_tx, paste0(outdir, sprintf("/results/truthtable_tx.csv")), sep="\t", row.names = F)



















