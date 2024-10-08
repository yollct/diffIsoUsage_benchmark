library(tidyverse)
library(ggplot2)
library(UpSetR)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(GGally)
library(fmsb)
library(gridExtra)
source("/nfs/home/students/chit/is_benchmark/plot_functions/functions.R")

path <- "/nfs/home/students/chit/is_benchmark"
outdir <- "/nfs/scratch/chit/simulated_real/%s_50_%s_%s%s"
noises <- c("", "_0.1","_0.5")

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

                        df_g <- df_g[grepl("ENSG", df_g$feature_id),]

                        outputpr <- cal_pre_re(df_g, truthfiles_g, split=split)
                        
                        this_pr <- pivot_output(outputpr, split=split)
                        this_pr$quant_tool <- quant
                        this_pr$rep <- paste(as.character(rep), 'replicates')
                        this_pr$seq <- seq
                        this_pr$sim <- R
                        this_pr$con <- paste0(quant, "_", paste(as.character(rep), 'replicates'))
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

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )


seq_label <- as_labeller(c('pair'="Paired-end", 'single'='Single-end', '2'='2', '3'='3','4'='4','5'='5', '2-4'='2-4','5-9'='5-9','>9'='>9','DTE'='S1','DTU'='S2','IS'='S3','0-0.1'='0-0.1', '0.1-0.2'='0.1-0.2','0.3-0.4'='0.3-0.4', '0.4-0.5'='0.4-0.5', '0.5-0.6'='0.5-0.6', '0.6-0.8'='0.6-0.8', '>0.8'='>0.8', '8 replicates'='8 replicates', '4 replicates'='4 replicates', '0'='0','0.5'='0.5', '0.1'='0.1', "kal" = "Kallisto", "rsem" = "RSEM","salmon" = "Salmon", 'kal_8 replicates'='Kallisto: 8 replicates', 'kal_4 replicates'='Kallisto: 4 replicates', 'rsem_8 replicates'='RSEM: 8 replicates', 'rsem_4 replicates'='RSEM: 4 replicates', 'salmon_8 replicates'='Salmon: 8 replicates', 'salmon_4 replicates'='Salmon: 4 replicates'))

alldf <- get_overall(split=NULL)
overalldf1 <- alldf %>% dplyr::filter(sim==whichrep&seq==seqtype) 

f1score <- overalldf1 %>% dplyr::filter(rep=="8 replicates"&noise==0.5)
write.table(f1score, "/nfs/scratch/chit/GSE222260/analysis/f1scores.csv", sep="\t")


plot_list <- c()
plot_name <- c()
n<-1
for (rps in c("4 replicates", "8 replicates")){
    for (spl in sort(unique(overalldf1$quant_tool))){
        data <- overalldf1 %>% filter(quant_tool==spl&rep==rps) %>% select(noise, f1, tool) %>% pivot_wider(names_from = tool, values_from = f1, id_cols = noise) %>% as.data.frame()
        row.names(data) <- paste0("bg-",data$noise)
        data <- rbind(rep(0.5,ncol(data)) , rep(0,ncol(data)) , data)
        data <- data %>% select(-noise)

        plot_list[[n]] <- data
        plot_name[[n]] <- paste0(rps,"_",spl)
        n<-n+1
 }
}

png(sprintf("./noise_fig/alltools_allf1radar_%s_%s_%s.png", seqtype, whichrep, thisrep), width=4000, height=2000, res=300)
nrows<-2
ncols<-3
par(mfrow = c(nrows, ncols), mar = c(0.5,0.5,0.5,0.5), oma=c(4,4,4,4)) # 4 rows, 2 columns

for (i in seq_along(plot_list)) {
radarchart(plot_list[[i]], axistype=1 , 
        #custom polygon
        pcol=colors_border ,  plwd=4 , plty=1,
        #custom the grid
        cglcol="grey", cglty=1, axislabcol="#5E5E5E", caxislabels=seq(0,0.6,0.15), cglwd=0.8,
        #custom labels
        vlcex=1.5, calcex=1
        
)  
}
row_labels <- c("4 replicates", "8 replicates") # Adjust as needed
for (i in 1:nrows) {
  mtext(row_labels[i], side = 2, at = (nrows - i + 0.5)/nrows, srt=90, outer = TRUE, line = 1, cex=1.8)
}

# Add column labels
col_labels <- c("Kallisto", "RSEM", "Salmon")  # Adjust as needed
for (i in 1:ncols) {
  mtext(col_labels[i], side = 3, at = (i - 0.5)/ncols, las = 1, outer = TRUE, line = 1, cex=1.8)
}

# 
dev.off()
