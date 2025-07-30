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
seqtype <- 'single'
whichrep <- 'r1'

thresholds <- data.frame(tools=c("iso_ktsp","drimseq","dexseq","dturtle","seqGSEA","cuffdiff","junctionseq","saturn","drimseq_stageR", "dexseq_stageR","saturn_stageR", "DSGseq", "nbsplice", "LimmaDS", "edgeR"),
                        thres=c(0.8,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05, 2, 0.05, 0.05, 0.05))
row.names(thresholds) <- thresholds$tools

get_table <- function(split){
    overalldf <- do.call(rbind, lapply(whichrep, function(R){
            do.call(rbind, lapply(seqtype, function(seq){
            do.call(rbind, lapply(thisrep, function(rep){
                do.call(rbind, lapply(noises, function(noise){
                    print(paste0(sprintf(outdir, seq, rep, R, noise)))
                    truthfiles_g <- read.csv(paste0(sprintf(outdir, seq, rep, R, noise), "/results/truthtable_gene.csv"), sep="\t") ### change to result
                    
                    
                    if (is.null(split)){
                        print("not consider DTE")
                        truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                    } else {
                        if (split!="events"){
                        print("not consider DTE")
                        truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
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

alldf <- get_table(split=NULL)


dsg_cut <- c()
isoktsp_cut <- c()
lapply(whichrep, function(R){
            lapply(seqtype, function(seq){
            lapply(thisrep, function(rep){
                lapply(noises, function(noise){
                    print(paste0(sprintf(outdir, seq, rep, R, noise)))
                    truthfiles_g <- read.csv(paste0(sprintf(outdir, seq, rep, R, noise), "/results/truthtable_gene.csv"), sep="\t") ### change to result
                    
                    
                    diffiso <- read.csv(paste0(sprintf(outdir, seq, rep, R, noise), "/results/isoinfo.csv"), sep="\t")
                    isorat <- diffiso %>% group_by(gene_id) %>% summarise(diffiso=max(diffiso))
                    truthfiles_g$diffiso <- isorat$diffiso
                    truthfiles_g <- truthfiles_g %>% dplyr::mutate(diffiso=ifelse(diffiso>0.8, ">0.8", ifelse(diffiso>=0.6, "0.6-0.8", ifelse(diffiso>=0.5, "0.5-0.6", ifelse(diffiso>=0.4, "0.4-0.5", ifelse(diffiso>=0.3, "0.3-0.4", ifelse(diffiso>=0.2, "0.2-0.3", ifelse(diffiso>=0.1, "0.1-0.2", "0-0.1")))))))) 
                    print(truthfiles_g %>% nrow)
                    # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    # truthfiles_tx <- truthfiles_tx %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    
                    lapply(c("salmon", "kal", "rsem"), function(quant){
                        df_g <- read.csv(paste0(sprintf(outdir, seq,rep, R, noise), "/results/", sprintf("%s_res_gene_N_T.txt", quant)), sep="\t")
                        thisthres <- thresholds[colnames(df_g)[2:ncol(df_g)],]$thres
                        df_g <- df_g %>% group_by(feature_id) %>% summarise(across(colnames(select(df_g, -feature_id)), .fns = min), .groups = "keep")
                        
                        
                        dsg_cut <<- c(dsg_cut, df_g$DSGseq)
                        isoktsp_cut <<- c(isoktsp_cut, df_g$iso_ktsp)
                        
                    })})
                })})
})

dsg_cut <- dsg_cut[!is.na(dsg_cut)]
dsg_cut <- dsg_cut[dsg_cut<5]
ggplot(data.frame(x=isoktsp_cut), aes(x))+
    geom_histogram(bins=10)


sal <- read.csv("/nfs/scratch/chit/simulated_real/pair_50_4_r1/results/salmon_res_gene_N_T.txt", sep="\t")

png("dsgdist_sim.png")
ggplot(data.frame(x=sal$DSGseq[sal$DSGseq]), aes(x))+
    geom_histogram(bins=30)+
    labs(x="NB stat distribution from DSGseq")+
    theme_light()+
    geom_vline(xintercept=2)+theme(axis.title = element_text(size=15))
dev.off()

png("isoktsp_sim.png")
ggplot(data.frame(x=sal$iso_ktsp), aes(x))+
    geom_histogram()+
    labs(x="Probability distribution from iso-KTSP")+
    theme_light()+
    geom_vline(xintercept=0.8)+
    theme(axis.title = element_text(size=15))
dev.off()

salr <- read.csv("/nfs/scratch/chit/GSE222260/results/salmon_res_gene_N_T.txt", sep="\t")
png("dsgdist.png")
ggplot(data.frame(x=salr$DSGseq[salr$DSGseq]), aes(x))+
    geom_histogram()+
    labs(x="NB stat distribution from DSGseq")+
    theme_light()+
    geom_vline(xintercept=5)+theme(axis.title = element_text(size=15))
dev.off()


#plot a false negative for a gene 
outdir <- "/nfs/scratch/chit/new_simulations/pair_50_8_r1_0"
gene <- "ENSG00000000005"

truthfile_tx <- read.csv(paste0(outdir, "/results/truthtable_tx.csv"), sep="\t")
truthfile_g <- read.csv(paste0(outdir, "/results/truthtable_gene.csv"), sep="\t")
salc <- read.csv("/nfs/scratch/chit/new_simulations/pair_50_8_r1_0/results/salmon_count.csv", sep=",")
salr <- read.csv(paste0(outdir,"/results/rsem_res_gene_N_T.txt"), sep="\t")
salc$feature_id <- lapply(salc$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist()
txid <- truthfile_tx[truthfile_tx$gene_id==gene,]$feature_id
salr[salr$feature_id %in% txid,]

truthfile_tx$feature_id %>% length
sum(salc$feature_id %in% truthfile_tx$feature_id)
sum(truthfile_tx$feature_id %in% salc$feature_id)
truthfile_tx %>% nrow
sum(truthfile_tx[truthfile_tx$status==1,]$feature_id %in% salc$feature_id)

truthfile_tx[truthfile_tx$status==1,]$feature_id %>% length
sum(truthfile_g$feature_id %in% salr$feature_id)
