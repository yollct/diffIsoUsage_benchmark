library(tidyverse)
library(ggplot2)
library(UpSetR)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(GGally)
source("/nfs/home/students/chit/is_benchmark/plot_functions/functions.R")

path <- "/nfs/home/students/chit/is_benchmark"
outdir <- "/nfs/scratch/chit/simulated_real/pair_%s_%s_%s%s"
dep <- c("50", "100")
noises <- c("_0.5")

thisrep <- c('4')
seqtype <- 'pair'
whichrep <- 'r1'

thresholds <- data.frame(tools=c("iso_ktsp","drimseq","dexseq","dturtle","seqGSEA","cuffdiff","junctionseq","saturn","drimseq_stageR", "dexseq_stageR","saturn_stageR", "DSGseq", "nbsplice", "LimmaDS", "edgeR"),
                        thres=c(0.5,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05, 5, 0.05, 0.05, 0.05))
row.names(thresholds) <- thresholds$tools

get_overall <- function(split){
    overalldf <- do.call(rbind, lapply(whichrep, function(R){
            do.call(rbind, lapply(dep, function(depth){
            do.call(rbind, lapply(thisrep, function(rep){
                do.call(rbind, lapply(noises, function(noise){
                    print(paste0(sprintf(outdir, depth, rep, R, noise)))
                    truthfiles_g <- read.csv(paste0(sprintf(outdir, depth, rep, R, noise), "/results/truthtable_gene.csv"), sep="\t") ### change to result
                    truthfiles_tx <- read.csv(paste0(sprintf(outdir, depth, rep, R, noise), "/results/truthtable_tx.csv"), sep="\t")
                    
                    
                    if (is.null(split)){
                        print("not consider DTE")
                        # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                    } else {
                        if (split!="events"){
                        print("not consider DTE")
                        # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                        }
                    }

                    diffiso <- read.csv(paste0(sprintf(outdir, depth, rep, R, noise), "/results/isoinfo.csv"), sep="\t")
                    isorat <- diffiso %>% group_by(gene_id) %>% summarise(diffiso=max(diffiso))
                    truthfiles_g$diffiso <- isorat$diffiso
                    truthfiles_g <- truthfiles_g %>% dplyr::mutate(diffiso=ifelse(diffiso>0.8, ">0.8", ifelse(diffiso>=0.6, "0.6-0.8", ifelse(diffiso>=0.5, "0.5-0.6", ifelse(diffiso>=0.4, "0.4-0.5", ifelse(diffiso>=0.3, "0.3-0.4", ifelse(diffiso>=0.2, "0.2-0.3", ifelse(diffiso>=0.1, "0.1-0.2", "0-0.1")))))))) 
                    print(truthfiles_g %>% nrow)
                    # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    # truthfiles_tx <- truthfiles_tx %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    
                    do.call(rbind, lapply(c("salmon", "kal", "rsem"), function(quant){
                        df_g <- read.csv(paste0(sprintf(outdir, depth,rep, R, noise), "/results/", sprintf("%s_res_gene_N_T.txt", quant)), sep="\t")
                        thisthres <- thresholds[colnames(df_g)[2:ncol(df_g)],]$thres
                        df_g <- df_g %>% group_by(feature_id) %>% summarise(across(colnames(select(df_g, -feature_id)), .fns = min), .groups = "keep")
                        

                        outputpr <- cal_pre_re(df_g, truthfiles_g, split=split)
                        
                        this_pr <- pivot_output(outputpr, split=split)
                        this_pr$quant_tool <- quant
                        this_pr$rep <- paste(as.character(rep), 'replicates')
                        this_pr$depth <- depth
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

seq_label <- as_labeller(c('pair'="Paired-end", 'single'='Single-end', '2'='2', '3'='3','4'='4','5'='5', '2-4'='2-4','5-9'='5-9','>9'='>9','DTE'='DTE','DTU'='DTU','IS'='IS','0-0.1'='0-0.1', '0.1-0.2'='0.1-0.2','0.3-0.4'='0.3-0.4', '0.4-0.5'='0.4-0.5', '0.5-0.6'='0.5-0.6', '0.6-0.8'='0.6-0.8', '>0.8'='>0.8', '8 replicates'='8 replicates', '4 replicates'='4 replicates', '0'='0','0.5'='0.5', '0.1'='0.1', "kal" = "Kallisto", "rsem" = "RSEM","salmon" = "Salmon", "50"="50", "100"="100"))

alldf <- get_overall(split=NULL)
overalldf1 <- alldf %>% dplyr::filter(sim==whichrep) 

overalldfkal <- overalldf1 %>% dplyr::filter(sim==whichrep&quant_tool=="sal")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("../noise_fig/depth_all_%s_%s_%s.png", "pair", whichrep, thisrep), width=900, height=600)
ggplot(overalldf1, aes(x=recall, y=precision, color=tool, shape=depth))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~quant_tool, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(colors="Tools", shape="Depth (Million)")
dev.off()
