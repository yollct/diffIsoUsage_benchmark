library(tidyverse)
library(ggplot2)
library(UpSetR)
source("/nfs/home/students/chit/is_benchmark/plot_functions/functions.R")

path <- "/nfs/home/students/chit/is_benchmark"
outdir <- "/nfs/scratch/chit/simulated_real/single4"
truthfiles_g <- list.files(paste0(outdir, "/results"), pattern="truthtable_gene.csv") ### change to result
truthfile_g <- read.csv(paste0(outdir, "/results/",truthfiles_g), sep="\t")
new_truthfile_g <- truthfile_g %>% select(-dispGeneEst) %>% unique()

truthfiles_tx <- list.files(paste0(outdir, "/results"), pattern="truthtable_tx.csv") ### change to result
truthfile_tx <- read.csv(paste0(outdir, "/results/",truthfiles_tx), sep="\t")

##
#truthfile_g$niso <- truthfile_tx %>% group_by(feature_id) %>% summarise(niso=n()) %>% dplyr::select(niso)
#truthfile_g <- truthfile_g %>% dplyr::mutate(gene_group=ifelse(niso<=1, 0, ifelse(niso>10, ">10", ifelse(niso>=6, "6-9", ifelse(niso>=3, "3-5", 2)))))

salmonfile <- list.files(paste0(outdir, "/results"), pattern="salmon_res_gene")
salmondf_g <- read.csv(paste0(outdir, "/results/", salmonfile), sep="\t")
salmondf_g[is.na(salmondf_g)] <- 1
salmondf_g <- salmondf_g %>% unique
salmondf_g <- salmondf_g[!grepl(".y",colnames(salmondf_g))]
colnames(salmondf_g) <- lapply(colnames(salmondf_g), function(x){gsub(".x", "", x)})

salmonfile <- list.files(paste0(outdir, "/results"), pattern="salmon_res_tx")
salmondf_tx <- read.csv(paste0(outdir, "/results/", kalfile), sep="\t")
salmondf_tx[is.na(salmondf_tx)] <- 1
salmon_tx$feature_id <- lapply(salmon_tx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist

kalfile <- list.files(paste0(outdir, "/results"), pattern="kal_res_gene")
kal_g <- read.csv(paste0(outdir, "/results/", kalfile), sep="\t")
kal_g[is.na(kal_g)] <- 1
kal_g <- kal_g %>% unique
kal_g <- kal_g[!grepl(".y",colnames(kal_g))]
colnames(kal_g) <- lapply(colnames(kal_g), function(x){gsub(".x", "", x)})

kalfile <- list.files(paste0(outdir, "/results"), pattern="kal_res_tx")
kal_tx <- read.csv(paste0(outdir, "/results/", kalfile), sep="\t")
kal_tx[is.na(kal_tx)] <- 1
kal_tx$feature_id <- lapply(kal_tx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist

rsemfile <- list.files(paste0(outdir, "/results"), pattern="rsem_res_gene")
rsemdf_g <- read.csv(paste0(outdir, "/results/", rsemfile), sep="\t")
rsemdf_g[is.na(rsemdf_g)] <- 1
rsemdf_g <- rsemdf_g %>% unique
rsemdf_g <- rsemdf_g[!grepl(".y", colnames(rsemdf_g))]

rsemfile <- list.files(paste0(outdir, "/results"), pattern="rsem_res_tx")
rsemdf_tx <- read.csv(paste0(outdir, "/results/", rsemfile), sep="\t")
rsemdf_tx[is.na(rsemdf_tx)] <- 1
####

sthresholds <- c(0.75, 0.05, 0.05,0.05, 0.05, 0.05,0.75,0.05,0.05)
kthresholds <- c(0.75,0.05,0.05,0.05,0.05,0.05, 0.75)

View(truthfile_tx)
outputpr<-cal_pre_re(kal_tx, truthfile_tx, kthresholds)
kal_pr <- pivot_output(outputpr)

png(paste0(outdir, "/results/kal_tx_pr_overall.png"))
ggplot(kal_pr, aes(x=precision, y=recall, color=tool, shape=tool, fill=tool))+
    geom_point(size=5, stroke=1.5)+
    scale_shape_manual(values=1:nlevels(kal_pr$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=20))
dev.off()
###
