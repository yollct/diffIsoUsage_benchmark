library(tidyverse)
library(ggplot2)
library(UpSetR)
source("/nfs/home/students/chit/is_benchmark/plot_functions/functions.R")

path <- "/nfs/home/students/chit/is_benchmark"
outdir <- "/nfs/scratch/chit/simulated_real/single_50_8_r1"

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
####

sthresholds <- thresholds[colnames(salmondf_g)[2:ncol(salmondf_g)],]$thres
print(sthresholds)
kthresholds <- thresholds[colnames(kal_g)[2:ncol(kal_g)],]$thres
print(kthresholds)
rthresholds <- thresholds[colnames(rsemdf_g)[2:ncol(rsemdf_g)],]$thres

res <- plot_upset(salmondf_g, sthresholds)
#salmondf_g <- salmondf_g %>% select(-iso_ktsp.x, -iso_ktsp.y)
res <- as.data.frame(res) 
res$feature_id <- salmondf_g$feature_id
colnames(res) <- lapply(colnames(res), function(x){strsplit(x, "[.]")[[1]][1]})
res <- res[,colnames(salmondf_g)]
png(paste0(outdir, "/results/upset.png"))
upset(res, nsets=7)
dev.off()

#####
# met <- salmondf_g$drimseq
# names(met) <- salmondf_g$feature_id
# names(met)[!new_truthfile_g[new_truthfile_g$status==1,]$feature_id %in% names(met)[met<0.05]]

outputpr<-cal_pre_re(salmondf_g, new_truthfile_g, split=NULL)

salmon_pr <- pivot_output(outputpr, split=NULL)
png(paste0(outdir, "/results/pr_overall.png"))
ggplot(salmon_pr, aes(x=precision, y=recall, color=tool, shape=tool, fill=tool))+
    geom_point(size=5, stroke=1.5)+
    scale_shape_manual(values=1:nlevels(salmon_pr$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=20))
dev.off()
###

outputpr<-cal_pre_re(salmondf_g, new_truthfile_g, split="gene_group")
salmon_pr_gg <- pivot_output(outputpr, split="gene_group")

png(paste0(outdir, "/results/pr_niso.png"), width=700)
ggplot(salmon_pr_gg, aes(x=precision, y=recall, color=tool, shape=tool))+
    geom_point(size=5, stroke=1.5) + facet_wrap(.~splits)+ggtitle("Grouped by number of switching isoforms")+
    scale_shape_manual(values=1:nlevels(salmon_pr_gg$tool))+
    theme_bw()+theme(axis.text=element_text(size=12), axis.title = element_text(size=20),strip.text.x = element_text(size = 20))
dev.off()
    
####

outputpr<-cal_pre_re(salmondf_g, new_truthfile_g, split="fc_group")
salmon_pr_fc <- pivot_output(outputpr, split="fc_group")

png(paste0(outdir, "/results/pr_foldchange.png"), width=700)
ggplot(salmon_pr_fc, aes(x=precision, y=recall, color=tool, shape=tool))+
    geom_point(size=5, stroke=1.5) + facet_wrap(.~splits)+ggtitle("Grouped by fold change")+
    scale_shape_manual(values=1:nlevels(salmon_pr_fc$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=12), axis.title = element_text(size=20))
dev.off()

#####
outputpr<-cal_pre_re(salmondf_g, new_truthfile_g, split="events")
salmon_pr_ev <- pivot_output(outputpr, split="events")

png(paste0(outdir, "/results/pr_events.png"),width=700)
ggplot(salmon_pr_ev, aes(x=precision, y=recall, color=tool, shape=tool))+
    geom_point(size=5, stroke=1.5) + facet_wrap(.~splits)+ggtitle("Grouped by fold change")+
    scale_shape_manual(values=1:nlevels(salmon_pr_ev$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=12), axis.title = element_text(size=20))
dev.off()


#################################### compare rsem with salmon
outputpr<-cal_pre_re(rsemdf_g, new_truthfile_g)
rsem_pr <- pivot_output(outputpr)

png(paste0(outdir, "/results/rsem_pr_overall.png"))
ggplot(rsem_pr, aes(x=precision, y=recall, color=tool, shape=tool, fill=tool))+
    geom_point(size=5, stroke=1.5)+
    scale_shape_manual(values=1:nlevels(rsem_pr$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=20))
dev.off()
###

outputpr<-cal_pre_re(rsemdf_g, new_truthfile_g, split="gene_group")
rsem_pr_gg <- pivot_output(outputpr, split="gene_group")

png(paste0(outdir, "/results/rsem_pr_niso.png"))
ggplot(rsem_pr_gg, aes(x=precision, y=recall, color=tool, shape=tool))+
    geom_point(size=5, stroke=1.5) + facet_wrap(.~splits)+ggtitle("Grouped by number of switching isoforms")+
    scale_shape_manual(values=1:nlevels(rsem_pr_gg$tool))+
    theme_bw()+theme(axis.text=element_text(size=12), axis.title = element_text(size=20))
dev.off()
    
####

outputpr<-cal_pre_re(rsemdf_g, new_truthfile_g, split="fc_group")
rsem_pr_fc <- pivot_output(outputpr, split="fc_group")

png(paste0(outdir, "/results/rsem_pr_foldchange.png"))
ggplot(rsem_pr_fc, aes(x=precision, y=recall, color=tool, shape=tool))+
    geom_point(size=5, stroke=1.5) + facet_wrap(.~splits)+ggtitle("Grouped by fold change")+
    scale_shape_manual(values=1:nlevels(rsem_pr_fc$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=12), axis.title = element_text(size=20))
dev.off()

#####
outputpr<-cal_pre_re(rsemdf_g, new_truthfile_g, split="events")
rsem_pr_ev <- pivot_output(outputpr, split="events")

png(paste0(outdir, "/results/rsem_pr_events.png"))
ggplot(rsem_pr_ev, aes(x=precision, y=recall, color=tool, shape=tool))+
    geom_point(size=5, stroke=1.5) + facet_wrap(.~splits)+ggtitle("Grouped by fold change")+
    scale_shape_manual(values=1:nlevels(rsem_pr_ev$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=12), axis.title = element_text(size=20))
dev.off()

##############
#kallisto

outputpr<-cal_pre_re(kal_g, new_truthfile_g)
kal_pr <- pivot_output(outputpr)

png(paste0(outdir, "/results/kal_pr_overall.png"))
ggplot(kal_pr, aes(x=precision, y=recall, color=tool, shape=tool, fill=tool))+
    geom_point(size=5, stroke=1.5)+
    scale_shape_manual(values=1:nlevels(kal_pr$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=15), axis.title = element_text(size=20))
dev.off()
###

outputpr<-cal_pre_re(kal_g, new_truthfile_g, split="gene_group")
kal_pr_gg <- pivot_output(outputpr, split="gene_group")

png(paste0(outdir, "/results/kal_pr_niso.png"), width=700)
ggplot(kal_pr_gg, aes(x=precision, y=recall, color=tool, shape=tool))+
    geom_point(size=5, stroke=1.5) + facet_wrap(.~splits)+ggtitle("Grouped by number of switching isoforms")+
    scale_shape_manual(values=1:nlevels(kal_pr_gg$tool))+
    theme_bw()+theme(axis.text=element_text(size=12), axis.title = element_text(size=20))
dev.off()
    
####

outputpr<-cal_pre_re(kal_g, new_truthfile_g, split="fc_group")
kal_pr_fc <- pivot_output(outputpr, split="fc_group")

png(paste0(outdir, "/results/kal_pr_foldchange.png"), width=700)
ggplot(kal_pr_fc, aes(x=precision, y=recall, color=tool, shape=tool))+
    geom_point(size=5, stroke=1.5) + facet_wrap(.~splits)+ggtitle("Grouped by fold change")+
    scale_shape_manual(values=1:nlevels(kal_pr_fc$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=12), axis.title = element_text(size=20))
dev.off()

#####
outputpr<-cal_pre_re(kal_g, new_truthfile_g, split="events")
kal_pr_ev <- pivot_output(outputpr, split="events")

png(paste0(outdir, "/results/kal_pr_events.png"), width=700)
ggplot(kal_pr_ev, aes(x=precision, y=recall, color=tool, shape=tool))+
    geom_point(size=5, stroke=1.5) + facet_wrap(.~splits)+ggtitle("Grouped by fold change")+
    scale_shape_manual(values=1:nlevels(kal_pr_ev$tool))+
    theme_bw()+
    theme(axis.text=element_text(size=12), axis.title = element_text(size=20))
dev.off()

png(paste0(outdir, "/results/kal_sal_pr.png"),width=700)
compare_tools(kal_pr, salmon_pr, rsem_pr)
dev.off()



















