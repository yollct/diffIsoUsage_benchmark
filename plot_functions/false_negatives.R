library(tidyverse)
library(ggplot2)
library(UpSetR)
library(rtracklayer)
source("/nfs/home/students/chit/is_benchmark/plot_functions/functions.R")

path <- "/nfs/home/students/chit/is_benchmark"
outdir <- "/nfs/scratch/chit/simulated_real/pair_50_8_r1/"

count <- read.csv(sprintf("%s/results/kal_count.csv", outdir), sep="\t")
#####
truthfiles_g <- list.files(paste0(outdir, "/results"), pattern="truthtable_gene.csv") ### change to result
truthfile_g <- read.csv(paste0(outdir, "/results/",truthfiles_g), sep="\t")
new_truthfile_g <- truthfile_g #%>% select(-dispGeneEst) %>% unique()

truthfiles_tx <- list.files(paste0(outdir, "/results"), pattern="truthtable_tx.csv") ### change to result
truthfile_tx <- read.csv(paste0(outdir, "/results/",truthfiles_tx), sep="\t")

salmonfile <- list.files(paste0(outdir, "/results"), pattern="salmon_res_gene")
salmondf_g <- read.csv(paste0(outdir, "/results/", salmonfile[1]), sep="\t")
salmondf_g[is.na(salmondf_g)] <- 1
salmondf_g <- salmondf_g %>% unique
# salmondf_g <- salmondf_g %>% select(-junctionseq)
# write.table(salmondf_g, paste0(outdir, "/results/", salmonfile[1]), sep="\t")
salmondf_g <- salmondf_g[!grepl(".y",colnames(salmondf_g))]
#colnames(salmondf_g) <- lapply(colnames(salmondf_g), function(x){gsub(".y", "", x)})

salmonfile <- list.files(paste0(outdir, "/results"), pattern="salmon_res_tx")
salmondf_tx <- read.csv(paste0(outdir, "/results/", salmonfile[1]), sep="\t")
salmondf_tx[is.na(salmondf_tx)] <- 1

kalfile <- list.files(paste0(outdir, "/results"), pattern="kal_res_gene")
kal_g <- read.csv(paste0(outdir, "/results/", kalfile[1]), sep="\t")
kal_g[is.na(kal_g)] <- 1
kal_g <- kal_g %>% unique

salmondf_g %>% dim
# kal_g <- kal_g[!grepl(".x",colnames(kal_g))]
colnames(kal_g) <- lapply(colnames(kal_g), function(x){gsub(".y$", "", x)})

kalfile <- list.files(paste0(outdir, "/results"), pattern="kal_res_tx")
kal_tx <- read.csv(paste0(outdir, "/results/", kalfile[1]), sep="\t")
kal_tx[is.na(kal_tx)] <- 1

rsemfile <- list.files(paste0(outdir, "/results"), pattern="rsem_res_gene")
rsemdf_g <- read.csv(paste0(outdir, "/results/", rsemfile[1]), sep="\t")
rsemdf_g[is.na(rsemdf_g)] <- 1
rsemdf_g <- rsemdf_g %>% unique
rsemdf_g <- rsemdf_g[!grepl(".y", colnames(rsemdf_g))]

rsemfile <- list.files(paste0(outdir, "/results"), pattern="rsem_res_tx")
rsemdf_tx <- read.csv(paste0(outdir, "/results/", rsemfile[1]), sep="\t")
rsemdf_tx[is.na(rsemdf_tx)] <- 1
####

gtf <- rtracklayer::import("/nfs/data/references/ensembl98_GRCh38/Homo_sapiens.GRCh38.98.gtf")
gtfgene <- as.data.frame(gtf) %>% dplyr::filter(type=="gene")

negs <- list_false_negative(salmondf_g, truthfile_g, "dexseq")
negs[1:10]

unfound <- negs[!negs %in% count$gene_id]
neggene <- gtfgene %>% filter(gene_id %in% unfound)

ggplot(neggene, aes(x=gene_biotype)) +
    geom_histogram(stat="count")

ggplot(neggene, aes(x=width)) +
    geom_histogram()

found <- negs[negs %in% count$gene_id]
neggene <- gtfgene %>% filter(gene_id %in% found)

posgene <- truthfile_g
truthfile_g %>% head
kal_g$feature_id[kal_g$cuffdiff>0.05] %in% truthfile_g$feature_id[truthfile_g$status==0]
