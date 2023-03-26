library(tidyverse)
library(ggplot2)
library(UpSetR)
source("/nfs/home/students/chit/is_benchmark/plot_functions/functions.R")

path <- "/nfs/home/students/chit/is_benchmark"
outdir <- "/nfs/scratch/chit/simulated_real/pair_50_4_r3/"

#####
truthfiles_g <- list.files(paste0(outdir, "/results"), pattern="truthtable_gene.csv") ### change to result
truthfile_g <- read.csv(paste0(outdir, "/results/",truthfiles_g), sep="\t")
new_truthfile_g <- truthfile_g #%>% select(-dispGeneEst) %>% unique()

truthfiles_tx <- list.files(paste0(outdir, "/results"), pattern="truthtable_tx.csv") ### change to result
truthfile_tx <- read.csv(paste0(outdir, "/results/",truthfiles_tx), sep="\t")

##
#truthfile_g$niso <- truthfile_tx %>% group_by(feature_id) %>% summarise(niso=n()) %>% dplyr::select(niso)
#truthfile_g <- truthfile_g %>% dplyr::mutate(gene_group=ifelse(niso<=1, 0, ifelse(niso>10, ">10", ifelse(niso>=6, "6-9", ifelse(niso>=3, "3-5", 2)))))

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

######
sal_stage <- read.csv(paste0(outdir, "/results/stager_salmon_res_gene_N_T.txt"), sep="\t")
sal_stage[is.na(sal_stage)] <- 1

kal_stage <- read.csv(paste0(outdir, "/results/stager_kal_res_gene_N_T.txt"), sep="\t")
kal_stage[is.na(kal_stage)] <- 1

rsem_stage <- read.csv(paste0(outdir, "/results/stager_rsem_res_gene_N_T.txt"), sep="\t")
rsem_stage[is.na(rsem_stage)] <- 1

sthresholds <- thresholds[colnames(salmondf_g)[2:ncol(salmondf_g)],]$thres
print(sthresholds)
kthresholds <- thresholds[colnames(kal_g)[2:ncol(kal_g)],]$thres
print(kthresholds)
rthresholds <- thresholds[colnames(rsemdf_g)[2:ncol(rsemdf_g)],]$thres

s_stage_thres <- thresholds[colnames(sal_stage)[2:ncol(sal_stage)],]$thres
k_stage_thres <- thresholds[colnames(kal_stage)[2:ncol(kal_stage)],]$thres
r_stage_thres <- thresholds[colnames(rsem_stage)[2:ncol(rsem_stage)],]$thres

res <- plot_upset(sal_stage, s_stage_thres)
#salmondf_g <- salmondf_g %>% select(-iso_ktsp.x, -iso_ktsp.y)
res <- as.data.frame(res) 
res$feature_id <- salmondf_g$feature_id
colnames(res) <- lapply(colnames(res), function(x){strsplit(x, "[.]")[[1]][1]})
res <- res[,colnames(salmondf_g)]
png(paste0(outdir, "/results/sstage_upset.png"))
upset(res, nsets=7)
dev.off()


outputpr<-cal_pre_re(sal_stage, new_truthfile_g, s_stage_thres, split=NULL)
ssalmon_pr <- pivot_output(outputpr, split=NULL)

png(paste0(outdir, "/results/stage_pr_overall.png"))
ggplot(salmon_pr, aes(x=precision, y=recall, color=tool, shape=tool, fill=tool))+
    geom_point(size=5, stroke=1.5)+
    scale_shape_manual(values=1:nlevels(salmon_pr$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=20))
dev.off()

outputpr<-cal_pre_re(rsem_stage, new_truthfile_g,  r_stage_thres)
rsem_pr <- pivot_output(outputpr)

png(paste0(outdir, "/results/stage_rsem_pr_overall.png"))
ggplot(rsem_pr, aes(x=precision, y=recall, color=tool, shape=tool, fill=tool))+
    geom_point(size=5, stroke=1.5)+
    scale_shape_manual(values=1:nlevels(rsem_pr$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=20))
dev.off()
###

outputpr<-cal_pre_re(kal_stage, new_truthfile_g,  k_stage_thres)
kal_pr <- pivot_output(outputpr)

png(paste0(outdir, "/results/stage_kal_pr_overall.png"))
ggplot(kal_pr, aes(x=precision, y=recall, color=tool, shape=tool, fill=tool))+
    geom_point(size=5, stroke=1.5)+
    scale_shape_manual(values=1:nlevels(kal_pr$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=20))
dev.off()
###


### diff genes after stageR 
before<-salmondf_g %>% filter(dexseq<0.05)
after<-sal_stage %>% filter(dexseq_stageR<0.05)
before %>% head

outputpr<-cal_pre_re(salmondf_g, new_truthfile_g, sthresholds, split=NULL)
salmon_pr <- pivot_output(outputpr, split=NULL)
