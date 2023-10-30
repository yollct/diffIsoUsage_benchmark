library(tidyverse)
library(ggplot2)
library(UpSetR)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(GGally)
source("/nfs/home/students/chit/is_benchmark/plot_functions/functions.R")

path <- "/nfs/home/students/chit/is_benchmark"
outdir <- "/nfs/scratch/chit/simulated_real/%s_50_%s_%s%s"
noises <- c("", "_0.1", "_0.5")

thisrep <- c('4','8')
seqtype <- 'pair'
whichrep <- 'r1'

thresholds <- data.frame(tools=c("iso_ktsp","drimseq","dexseq","dturtle","seqGSEA","cuffdiff","junctionseq","saturn","drimseq_stageR", "dexseq_stageR","saturn_stageR", "DSGseq", "nbsplice", "LimmaDS", "edgeR"),
                        thres=c(0.5,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05, 5, 0.05, 0.05, 0.05))
row.names(thresholds) <- thresholds$tools

get_overall <- function(split){
    overalldf <- do.call(rbind, lapply(whichrep, function(R){
            do.call(rbind, lapply(seqtype, function(seq){
            do.call(rbind, lapply(thisrep, function(rep){
                do.call(rbind, lapply(noises, function(noise){
                    print(paste0(sprintf(outdir, seq, rep, R, noise)))
                    truthfiles_g <- read.csv(paste0(sprintf(outdir, seq, rep, R, noise), "/results/truthtable_gene.csv"), sep="\t") ### change to result
                    truthfiles_tx <- read.csv(paste0(sprintf(outdir, seq, rep, R, noise), "/results/truthtable_tx.csv"), sep="\t")
                    
                    
                    if (is.null(split)){
                        print("not consider DTE")
                        # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                    } else {
                        if (split!="events"){
                        print("not consider DTE")
                        # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                        }
                    }

                    diffiso <- read.csv(paste0(sprintf(outdir, seq, rep, R, noise), "/results/isoinfo.csv"), sep="\t")
                    isorat <- diffiso %>% group_by(gene_id) %>% summarise(diffiso=max(diffiso))
                    truthfiles_g$diffiso <- isorat$diffiso
                    truthfiles_g <- truthfiles_g %>% dplyr::mutate(diffiso=ifelse(diffiso>0.8, ">0.8", ifelse(diffiso>=0.6, "0.6-0.8", ifelse(diffiso>=0.5, "0.5-0.6", ifelse(diffiso>=0.4, "0.4-0.5", ifelse(diffiso>=0.3, "0.3-0.4", ifelse(diffiso>=0.2, "0.2-0.3", ifelse(diffiso>=0.1, "0.1-0.2", "0-0.1")))))))) 
                    print(truthfiles_g %>% nrow)
                    # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    # truthfiles_tx <- truthfiles_tx %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    
                    do.call(rbind, lapply(c("salmon", "kal", "rsem"), function(quant){
                        df_g <- read.csv(paste0(sprintf(outdir, seq,rep, R, noise), "/results/", sprintf("%s_res_gene_N_T.txt", quant)), sep="\t")
                        thisthres <- thresholds[colnames(df_g)[2:ncol(df_g)],]$thres
                        df_g <- df_g %>% group_by(feature_id) %>% summarise(across(colnames(select(df_g, -feature_id)), .fns = min), .groups = "keep")
                        

                        outputpr <- cal_pre_re(df_g, truthfiles_g, split=split)
                        
                        this_pr <- pivot_output(outputpr, split=split)
                        this_pr$quant_tool <- quant
                        this_pr$rep <- paste(as.character(rep), 'replicates')
                        this_pr$seq <- seq
                        this_pr$sim <- R
                        this_pr$noise <- ifelse(noise=="", 0, ifelse(noise=="_0.1", 0.1, 0.5))
                        rm(df_g)
                        this_pr
                    }))
                }))
            }))
        }))
    }))
    overalldf$f1 <- lapply(overalldf$f1, function(x){as.numeric(x)}) %>% unlist
    overalldf
}

seq_label <- as_labeller(c('pair'="Paired-end", 'single'='Single-end', '2'='2', '3'='3','4'='4','5'='5', '2-4'='2-4','5-9'='5-9','>9'='>9','DTE'='DTE','DTU'='DTU','IS'='IS','0-0.1'='0-0.1', '0.1-0.2'='0.1-0.2','0.3-0.4'='0.3-0.4', '0.4-0.5'='0.4-0.5', '0.5-0.6'='0.5-0.6', '0.6-0.8'='0.6-0.8', '>0.8'='>0.8', '8 replicates'='8 replicates', '4 replicates'='4 replicates', '0'='0','0.5'='0.5', '0.1'='0.1', "kal" = "Kallisto", "rsem" = "RSEM","salmon" = "Salmon"))

alldf <- get_overall(split=NULL)
overalldf %>% head
overalldf1 <- alldf %>% dplyr::filter(sim==whichrep&seq==seqtype) 

f1score <- overalldf1 %>% dplyr::filter(rep=="8 replicates"&noise==0.5)
write.table(f1score, "/nfs/scratch/chit/GSE222260/analysis/f1scores.csv", sep="\t")

png(sprintf("./noise_fig/alltools_re_%s_%s_%s.png", seqtype, whichrep, thisrep), width=1000, height=600)
ggplot(overalldf1, aes(x=noise, y=recall, color=tool))+
    geom_point(size=7, alpha=0.8)+geom_path()+
    scale_color_manual(values=colMap)+
  facet_grid(rep~quant_tool, labeller=seq_label)+
    theme(axis.text=element_text(size=16), axis.text.x = element_text(vjust = 0, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Recall")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")+theme_light()
dev.off()

png(sprintf("./noise_fig/alltools_pr_%s_%s_%s.png", seqtype, whichrep, thisrep), width=1000, height=600)
ggplot(overalldf1, aes(x=noise, y=precision, color=tool))+
    geom_point(size=7, alpha=0.8)+geom_path()+
    scale_color_manual(values=colMap)+
  facet_grid(rep~quant_tool, labeller=seq_label)+
    theme(axis.text=element_text(size=16), axis.text.x = element_text(vjust = 0, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Precision")+
    labs(color="Tools", shape="Quantification tools")+theme_light()
dev.off()

png(sprintf("./noise_fig/alltools_f1_%s_%s_%s.png", seqtype, whichrep, thisrep), width=1000, height=600)
ggplot(overalldf1, aes(x=noise, y=f1, color=tool))+
    geom_point(size=7, alpha=0.8)+geom_path()+
    scale_color_manual(values=colMap)+
  facet_grid(rep~quant_tool, labeller=seq_label)+
    theme(axis.text=element_text(size=16), axis.text.x = element_text(vjust = 0, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Precision")+
    labs(color="Tools", shape="Quantification tools")+theme_light()
dev.off()

overalldfkal <- overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/alltools_all_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(overalldf1, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~quant_tool, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(colors="Tools", shape="Invariance prob.")
dev.off()

overalldfsal <- alldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=500, height=600)
ggplot(overalldfsal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~., labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Noise")
dev.off()

fc_overalldf <- get_overall(split="fc_group")
fc_overalldf1 <- fc_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype)
write.table(fc_overalldf, "./fc_overalldf_pair.csv", sep="\t", row.names=F)

png(sprintf("./noise_fig/fc_re_alltools_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(fc_overalldf1, aes(x=splits, y=recall, color=tool, shape=quant_tool, alpha=rep))+
    geom_path()+
    geom_point(size=7, alpha=0.8)+ 
   scale_color_manual(values=colMap)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Recall")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")+theme_light()
dev.off()

png(sprintf("./noise_fig/fc_alltools_pr_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(fc_overalldf1, aes(x=noise, y=precision, color=tool, shape=quant_tool))+
    geom_path()+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Precision")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")+theme_light()
dev.off()

fc_overalldfkal <- fc_overalldf1 <- fc_overalldf %>% dplyr::filter(quant_tool=="kal")
png(sprintf("./noise_fig/fc_alltools_0.5_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(fc_overalldfkal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=7, alpha=0.8)+theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Precision")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")+theme_light()
dev.off()

fc_overalldfkal <- fc_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/fc_alltools_kal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(fc_overalldfkal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()

fc_overalldfsal <- fc_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/fc_alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(fc_overalldfsal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()

ir_overalldf <- get_overall(split="diffiso")
ir_overalldf1 <- ir_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&splits!='0.6-0.8'&splits!='>0.8')


png(sprintf("./noise_fig/diffiso_re_alltools_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(ir_overalldf1, aes(x=noise, y=recall, color=tool, shape=quant_tool))+
    geom_path()+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Recall")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()


png(sprintf("./noise_fig/ir_alltools_pr_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(ir_overalldf1, aes(x=noise, y=precision, color=tool, shape=quant_tool))+
    geom_path()+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Precision")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()

ir_overalldfkal <- ir_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/ir_alltools_kal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(ir_overalldfkal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()

ir_overalldfsal <- ir_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/ir_alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(ir_overalldfsal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()

iso_overalldf <- get_overall(split="diffiso_group")
iso_overalldf1 <- iso_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype)

png(sprintf("./noise_fig/diffisogroup_re_alltools_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(iso_overalldf1, aes(x=noise, y=recall, color=tool, shape=quant_tool))+
    geom_path()+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Recall")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()


png(sprintf("./noise_fig/diffisogroup_alltools_pr_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(iso_overalldf1, aes(x=noise, y=precision, color=tool, shape=quant_tool))+
    geom_path()+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Precision")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()

iso_overalldfkal <- iso_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/diffiso_alltools_kal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(iso_overalldfkal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()

iso_overalldfsal <- iso_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/diffiso_alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(iso_overalldfsal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()

ev_overalldf <- get_overall(split="events")
ev_overalldf1 <- ev_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype)

png(sprintf("./noise_fig/ev_re_alltools_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(ev_overalldf1, aes(x=noise, y=recall, color=tool, shape=quant_tool))+
    geom_path()+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Recall")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()


png(sprintf("./noise_fig/ev_alltools_pr_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(ev_overalldf1, aes(x=noise, y=precision, color=tool, shape=quant_tool))+
    geom_path()+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Precision")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()

ev_overalldfkal <- ev_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/ev_alltools_kal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(ev_overalldfkal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()

ev_overalldfsal <- ev_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/ev_alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(ev_overalldfsal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()

gg_overalldf <- get_overall(split="gene_group")
gg_overalldf1 <- gg_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype)

png(sprintf("./noise_fig/genegroup_re_alltools_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(gg_overalldf1, aes(x=noise, y=recall, color=tool, shape=quant_tool))+
    geom_path()+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Recall")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()


png(sprintf("./noise_fig/genegroup_alltools_pr_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(gg_overalldf1, aes(x=noise, y=precision, color=tool, shape=quant_tool))+
    geom_path()+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Precision")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()


gg_overalldfkal <- gg_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/genegroup_alltools_kal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(gg_overalldfkal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()

gg_overalldfsal <- gg_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./noise_fig/genegroup_alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(gg_overalldfsal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()



# tru <- tru %>% dplyr::filter(fc_group==4)
# p_g <- tru$feature_id[tru$status==1]
# dex1 <- dex %>% dplyr::mutate(dexseq=replace_na(dexseq, 1))
# dex1 <- dex1 %>% group_by(feature_id) %>% summarise(across(colnames(select(dex, -feature_id)), .fns = min), .groups = "keep")

# d_g <- dex1$feature_id[dex1$dexseq<0.05]

# sal <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/rsem_count.csv")
# sum(is.na(dex$feature_id))
tru <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r2_0.1/results/truthtable_gene.csv", sep="\t")

diffiso <- read.csv("/nfs/scratch/chit/simulated_real/pair_50_4_r1/results/isoinfo.csv", sep="\t")
isorat <- diffiso %>% group_by(gene_id) %>% summarise(diffiso=max(diffiso))
tru$diffiso <- isorat$diffiso
tru <- tru %>% dplyr::mutate(diffiso=ifelse(diffiso>0.8, ">0.8", ifelse(diffiso>=0.6, "0.6-0.8", ifelse(diffiso>=0.5, "0.5-0.6", ifelse(diffiso>=0.4, "0.4-0.5", ifelse(diffiso>=0.3, "0.3-0.4", ifelse(diffiso>=0.2, "0.2-0.3", ifelse(diffiso>=0.1, "0.1-0.2", "0-0.1")))))))) 
tru <- tru %>% dplyr::filter(fc_group!="1")
png("./noise_fig/tru_stat.png")
ggplot(tru, aes(x=diffiso))+
    geom_bar()+
    facet_grid(gene_group~fc_group)
dev.off()
# sum(d_g %in% p_g)/length(p_g)
sal <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r2_0.1/results/salmon_res_gene_N_T.txt", sep="\t")

sal$drimseq[is.na(sal$drimseq)] <- 1
namet <- sal$feature_id[sal$drimseq<0.05]
nrow(sal)
fctru <- tru %>% filter(events=="IS")
fctru %>% head
fctru %>% nrow

#not_tru <- tru %>% filter(fc_group!="3") 

fil1 <- sal[sal$feature_id %in% tru$feature_id,]

detectedpos <- fil1$feature_id[fil1$drimseq <0.05]
boomet <- unique(detectedpos)
tru <- tru %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))

split_tru <- tru %>% filter(status==1)
split_tru %>% nrow
detectedpos[detectedpos %in% split_tru$feature_id] 

tru %>% head

fc_overalldf %>% filter(tool=="nbsplice" & quant_tool=="salmon")
