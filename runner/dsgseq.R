library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
path <- args[2]
con1 <- args[3]
con2 <- args[4]
print(con1)
print(con2)

dsg_res <- read_tsv(paste0(outdir, sprintf("/results/DSGseq_results.diff")))

resgene <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)))
kresgene <- read_tsv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)))
rsresgene <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)))

restx <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)))
krestx <- read_tsv(paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)))
rsrestx <- read_tsv(paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)))

dsg_res1 <- dsg_res %>% group_by(gene_name) %>% summarise(NB_stat=max(NB_stat))
dsg_g <- dsg_res1 %>% dplyr::select(gene_name, NB_stat) %>% dplyr::rename('feature_id'='gene_name', 'DSGseq'='NB_stat')
dsg_tx <- dsg_res %>% dplyr::select(ID, NB_stat) %>% dplyr::rename('feature_id'="ID", "DSGseq"="NB_stat")

if (!any(grepl("DSGSeq", colnames(resgene)))) { 
    resgene <- full_join(resgene, dsg_g, by="feature_id")
    restx <- full_join(restx, dsg_tx, by="feature_id")
}

if (!any(grepl("DSGSeq", colnames(kresgene)))) { 
    kresgene <- full_join(kresgene, dsg_g, by="feature_id")
    krestx <- full_join(krestx, dsg_tx, by="feature_id")
}

if (!any(grepl("DSGSeq", colnames(rsresgene)))) { 
    rsresgene <- full_join(rsresgene, dsg_g, by="feature_id")
    rsrestx <- full_join(rsrestx, dsg_tx, by="feature_id")
}

write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")

write.table(kresgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")

write.table(rsresgene, paste0(outdir, sprintf("/results/rsem_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")


write.table(restx, paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")

write.table(krestx, paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")

write.table(rsrestx, paste0(outdir, sprintf("/results/rsem_res_tx_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")


