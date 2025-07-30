library(tidyverse)
library(ggplot2)
library(UpSetR)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(GGally)
library(fmsb)
library(gridExtra)
source("/home/chit/norm_is/plot_functions/new_functions.R")

path <- "/home/chit/norm_is"
outdir <- "/home/chit/norm_chit/new_simulations/%s_50_%s_%s%s"
noises <- c("_0", "_0.1","_0.5", "_0.7", "_0.9")

thisrep <- c('4','8')
seqtype <- 'pair'
whichrep <- 'r1'

thresholds <- data.frame(tools=c("iso_ktsp","drimseq","dexseq","dturtle","seqGSEA","cuffdiff","junctionseq","saturn","drimseq_stageR", "dexseq_stageR","saturn_stageR", "DSGseq", "nbsplice", "LimmaDS", "edgeR", "BANDITS"),
                        thres=c(0.5,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05, 5, 0.05, 0.05, 0.05, 0.05))
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

                    # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    # truthfiles_tx <- truthfiles_tx %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    
                    do.call(rbind, lapply(c("salmon", "kal", "rsem"), function(quant){
                        df_g <- read.csv(paste0(sprintf(outdir, seq,rep, R, noise), "/results/", sprintf("%s_res_tx_N_T.txt", quant)), sep="\t")
                        print("read")
                        thisthres <- thresholds[colnames(df_g)[2:ncol(df_g)],]$thres
                        df_g <- df_g %>% group_by(feature_id) %>% summarise(across(colnames(select(df_g, -feature_id)), .fns = min), .groups = "keep")

                        df_g <- df_g[grepl("ENST", df_g$feature_id),]
                        this_tru <- truthfiles_tx[truthfiles_tx$feature_id %in% df_g$feature_id,]

                        
                        outputpr <- cal_pre_re(df_g, this_tru, split=split)
                        print("calculated")
                        this_pr <- pivot_output(outputpr, split=split)
                        this_pr$quant_tool <- quant
                        this_pr$rep <- paste(as.character(rep), 'replicates')
                        this_pr$seq <- seq
                        this_pr$sim <- R
                        this_pr$con <- paste0(quant, "_", paste(as.character(rep), 'replicates'))
                        this_pr$noise <- gsub("_", "", noise)
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



seq_label <- as_labeller(c('pair'="Paired-end", 'single'='Single-end', '2'='2', '3'='3','4'='4','5'='5', '2-4'='2-4','5-8'='5-8','>9'='>9','DTE'='S1','DTU'='S2','IS'='S3','0-0.1'='0-0.1', '0.1-0.2'='0.1-0.2','0.3-0.4'='0.3-0.4', '0.4-0.5'='0.4-0.5', '0.5-0.6'='0.5-0.6', '0.6-0.8'='0.6-0.8', '>0.8'='>0.8', '8 replicates'='8 replicates', '4 replicates'='4 replicates', '0'='0','0.5'='0.5', '0.1'='0.1', '0.7'='0.7', '0.9'='0.9',"kal" = "Kallisto", "rsem" = "RSEM","salmon" = "Salmon", 'kal_8 replicates'='Kallisto: 8 replicates', 'kal_4 replicates'='Kallisto: 4 replicates', 'rsem_8 replicates'='RSEM: 8 replicates', 'rsem_4 replicates'='RSEM: 4 replicates', 'salmon_8 replicates'='Salmon: 8 replicates', 'salmon_4 replicates'='Salmon: 4 replicates'))

alldf <- get_overall(split=NULL)
write.table(alldf, "/home/chit/norm_is/alldf.csv", sep="\t")

q()
overalldf1 <- alldf %>% dplyr::filter(sim==whichrep&seq==seqtype) 

f1score <- overalldf1 %>% dplyr::filter(rep=="8 replicates"&noise==0.5)
write.table(f1score, "/nfs/scratch/chit/GSE222260/analysis/f1scores.csv", sep="\t")

# png(sprintf("./noise_fig/alltools_re_%s_%s_%s.png", seqtype, whichrep, thisrep), width=1000, height=600)
# ggplot(overalldf1, aes(x=noise, y=recall, color=tool))+
#     geom_point(size=7, alpha=0.8)+geom_path()+
#     scale_color_manual(values=colMap)+
#   facet_grid(rep~quant_tool, labeller=seq_label)+
#     theme(axis.text=element_text(size=16), axis.text.x = element_text(vjust = 0, hjust = 1), axis.title = element_text(size=20),
#         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
#     xlab("Noise level")+
#     ylab("Recall")+
#     scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
#                               "salmon" = "Salmon"))+
#     labs(color="Tools", shape="Quantification tools")+theme_light()
# dev.off()

# png(sprintf("./noise_fig/alltools_pr_%s_%s_%s.png", seqtype, whichrep, thisrep), width=1000, height=600)
# ggplot(overalldf1, aes(x=noise, y=precision, color=tool))+
#     geom_point(size=7, alpha=0.8)+geom_path()+
#     scale_color_manual(values=colMap)+
#   facet_grid(rep~quant_tool, labeller=seq_label)+
#     theme(axis.text=element_text(size=16), axis.text.x = element_text(vjust = 0, hjust = 1), axis.title = element_text(size=20),
#         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
#     xlab("Noise level")+
#     ylab("Precision")+
#     labs(color="Tools", shape="Quantification tools")+theme_light()
# dev.off()

jpeg(sprintf("%s/tx_noise_fig/alltools_allf1_%s_%s.jpeg", path, seqtype, whichrep), res=300, width=3000, height=2000)
overalldf1 %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>%
ggplot(aes(x=tool, y=f1, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=f1, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~quant_tool, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
        ylab("F1 score")+xlab("Tool")

dev.off()

jpeg(sprintf("%s/tx_noise_fig/alltools_f1_%s_%s.jpeg", path,seqtype, whichrep), width=4000, height=2000)
ggplot(overalldf1, aes(x=noise, y=f1, color=tool))+
    geom_point(size=7, alpha=0.8)+geom_path()+
    scale_color_manual(values=colMap)+
  facet_grid(rep~quant_tool, labeller=seq_label)+
    theme(axis.text=element_text(size=16), axis.text.x = element_text(vjust = 0, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Background")+
    ylab("Precision")+
    labs(color="Tools", shape="Quantification tools")+theme_light()
dev.off()

overalldfkal <- overalldf1 %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
jpeg(sprintf("%s/tx_noise_fig/alltools_all_%s_%s.jpeg", path,seqtype, whichrep), res=300, width=4000, height=2000)
ggplot(overalldf1, aes(x=recall, y=as.numeric(fdr), color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=as.numeric(fdr), group=tool))+
    theme_light()+
    labs(colors="Tools", shape="Background",xlab="TPR", ylab="FDR") +
   scale_color_manual(values=colMap)+
   facet_grid(rep~quant_tool, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) 
dev.off()

overalldfsal <- alldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
jpeg(sprintf("%s/tx_noise_fig/alltools_sal_%s_%s.jpeg", path, seqtype, whichrep), res=300, width=4000, height=2000)
ggplot(overalldfsal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~., labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
         labs(xlab="Recall", ylab="Precision")
dev.off()

png(sprintf("./tx_noise_fig/alltools_%s_%s_%s_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
ggplot(overalldf1, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(con~., labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
         labs(xlab="Recall", ylab="Precision")
dev.off()


fc_overalldf <- get_overall(split="fc_group")

fc_overalldf1 <- fc_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype)
fc_overalldf1$fdr[is.na(fc_overalldf1$fdr)] <- 0
#radar plot


plot_list <- c()
plot_name <- c()
n<-1
for (rps in c("4 replicates", "8 replicates")){
    for (spl in sort(unique(fc_overalldf1$splits))){
        data <- fc_overalldf1 %>% filter(quant_tool=="salmon"&splits==spl&rep==rps) %>% select(noise, fdr, tool) %>% pivot_wider(names_from = tool, values_from = fdr, id_cols = noise) %>% as.data.frame()
        row.names(data) <- paste0("bg-",data$noise)
        data <- rbind(rep(1,ncol(data)) , rep(0,ncol(data)) , data)
        data <- data %>% select(-noise)
        data[is.na(data)] <- 0
        data <- lapply(data, function(x) if(is.factor(x) | is.character(x)) as.numeric(as.character(x)) else x)

        plot_list[[n]] <- as.data.frame(data)
        plot_name[[n]] <- paste0(rps,"_",spl)
        n<-n+1
 }
}

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9), "#008000" , "black" )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4),"#008000" , "black" )

png(sprintf("./tx_noise_fig/fc_allradar_%s_%s_%s.png", seqtype, whichrep, thisrep), width=4500, height=2000, res=300)
nrows<-2
ncols<-4
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
col_labels <- sort(unique(fc_overalldf1$splits))   # Adjust as needed
for (i in 1:ncols) {
  mtext(col_labels[i], side = 3, at = (i - 0.5)/ncols, las = 1, outer = TRUE, line = 1, cex=1.8)
}

# 
dev.off()

png("./tx_noise_fig/radar_legend.png", width = 2000, height = 2000, res=300)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
legend(x=0.5, y=0.5, legend = c("0","0.1","0.5"), title="Background", bty = "n", pch=20 , col=colors_border , text.col = "#494949", cex=2, pt.cex=4)
dev.off()

fc_overalldfkal <- fc_overalldf1 <- fc_overalldf %>% dplyr::filter(quant_tool=="kal")
png(sprintf("./tx_noise_fig/fc_alltools_0.5_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ggplot(fc_overalldfkal, aes(x=recall, y=precision, shape=as.character(noise)))+
    geom_point(size=7, alpha=0.8)+theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Noise level")+
    ylab("Precision")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools", shape="Background")+theme_light()
dev.off()

fc_overalldfkal <- fc_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal") %>% as.data.frame
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./tx_noise_fig/fc_alltools_kal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
fc_overalldfkal %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>%
ggplot(aes(x=tool, y=f1, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        labs(color="Tools", shape="Background")
dev.off()

fc_overalldfsal <- fc_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./tx_noise_fig/fc_alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
fc_overalldfsal %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>%
ggplot(aes(x=tool, y=f1, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   ylim(c(0,0.6))+
   theme(axis.text=element_text(size=16), axis.text.x = element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        labs(color="Tools", shape="Background.")
dev.off()

png(sprintf("./tx_noise_fig/fc_alltools_%s_%s_%s_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
fc_overalldf %>%
    mutate(tool=fct_reorder(tool, desc(f1)), fdr=as.numeric(fdr)) %>%
ggplot(aes(x=tool, y=fdr, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=f1, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(con~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
        ggtitle("FDR in different event scenarios in paired-end data")
dev.off()

png(sprintf("./tx_noise_fig/fc_alltools_%s_%s_%s_re_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
fc_overalldf %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>%
ggplot(aes(x=tool, y=recall, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=recall, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(con~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
        ggtitle("Recall in different event scenarios in paired-end data")
dev.off()


png(sprintf("./tx_noise_fig/fc_alltools_%s_%s_%s_pr_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
fc_overalldf %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>%
ggplot(aes(x=tool, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(con~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
        ggtitle("Precision in different fold change in paired-end data")
dev.off()

# ir_overalldf <- get_overall(split="diffiso")
# ir_overalldf1 <- ir_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&splits!='0.6-0.8'&splits!='>0.8')


# # png(sprintf("./noise_fig/diffiso_re_alltools_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
# # ggplot(ir_overalldf1, aes(x=noise, y=recall, color=tool, shape=quant_tool))+
# #     geom_path()+
# #     geom_point(size=7, alpha=0.8)+ theme_light()+
# #    scale_color_manual(values=colMap)+
# #    facet_grid(rep~splits, labeller=seq_label)+
# #    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
# #         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
# #     xlab("Noise level")+
# #     ylab("Recall")+
# #     scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
# #                               "salmon" = "Salmon"))+
# #     labs(color="Tools", shape="Quantification tools")
# # dev.off()


# # png(sprintf("./noise_fig/ir_alltools_pr_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
# # ggplot(ir_overalldf1, aes(x=noise, y=precision, color=tool, shape=quant_tool))+
# #     geom_path()+
# #     geom_point(size=7, alpha=0.8)+ theme_light()+
# #    scale_color_manual(values=colMap)+
# #    facet_grid(rep~splits, labeller=seq_label)+
# #    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
# #         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
# #     xlab("Noise level")+
# #     ylab("Precision")+
# #     scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
# #                               "salmon" = "Salmon"))+
# #     labs(color="Tools", shape="Quantification tools")
# # dev.off()

# ir_overalldfkal <- ir_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal")
# # dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
# png(sprintf("./tx_noise_fig/ir_alltools_kal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
# ggplot(ir_overalldfkal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
#     geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
#     theme_light()+
#    scale_color_manual(values=colMap)+
#    facet_grid(rep~splits, labeller=seq_label)+
#    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
#         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
#         labs(color="Tools", shape="Background")
# dev.off()

# ir_overalldfsal <- ir_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# # dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
# png(sprintf("./tx_noise_fig/ir_alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
# ggplot(ir_overalldfsal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
#     geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
#     theme_light()+
#    scale_color_manual(values=colMap)+
#    facet_grid(rep~splits, labeller=seq_label)+
#    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
#         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
#         labs(color="Tools", shape="Background")
# dev.off()

# png(sprintf("./tx_noise_fig/ir_alltools_%s_%s_%s_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
# ggplot(ir_overalldf, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
#     geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
#     theme_light()+
#    scale_color_manual(values=colMap)+
#    facet_grid(con~splits, labeller=seq_label)+
#    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
#         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
#         labs(color="Tools", shape="Background")
# dev.off()

# iso_overalldf <- get_overall(split="diffiso_group")
# iso_overalldf1 <- iso_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype)

# #radar plot
# library(fmsb)
# library(gridExtra)

# plot_list <- c()
# plot_name <- c()
# n<-1
# for (rps in c("4 replicates", "8 replicates")){
#     for (spl in sort(unique(iso_overalldf1$splits))){
#         data <- iso_overalldf1 %>% filter(quant_tool=="salmon"&splits==spl&rep==rps) %>% select(noise, f1, tool) %>% pivot_wider(names_from = tool, values_from = f1, id_cols = noise) %>% as.data.frame()
#         row.names(data) <- paste0("bg-",data$noise)
#         data <- rbind(rep(0.5,ncol(data)) , rep(0,ncol(data)) , data)
#         data <- data %>% select(-noise)

#         plot_list[[n]] <- data
#         plot_name[[n]] <- paste0(rps,"_",spl)
#         n<-n+1
#  }
# }

# colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9), "#008000" , "black" )
# colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4),"#008000" , "black" )

# png(sprintf("./tx_noise_fig/iso_allradar_%s_%s_%s.png", seqtype, whichrep, thisrep), width=4000, height=2000, res=300)
# nrows<-2
# ncols<-4
# par(mfrow = c(nrows, ncols), mar = c(0.5,0.5,0.5,0.5), oma=c(4,4,4,4)) # 4 rows, 2 columns

# for (i in seq_along(plot_list)) {
# radarchart(plot_list[[i]], axistype=1 , 
#         #custom polygon
#         pcol=colors_border ,  plwd=4 , plty=1,
#         #custom the grid
#         cglcol="grey", cglty=1, axislabcol="#5E5E5E", caxislabels=seq(0,0.5,0.1), cglwd=0.8,
#         #custom labels
#         vlcex=1.5, calcex=1.0
        
# )  
# }
# row_labels <- c("4 replicates", "8 replicates") # Adjust as needed
# for (i in 1:nrows) {
#   mtext(row_labels[i], side = 2, at = (nrows - i + 0.5)/nrows, srt=90, outer = TRUE, line = 1, cex=1.8)
# }

# # Add column labels
# col_labels <- sort(unique(iso_overalldf1$splits))   # Adjust as needed
# for (i in 1:ncols) {
#   mtext(col_labels[i], side = 3, at = (i - 0.5)/ncols, las = 1, outer = TRUE, line = 1, cex=1.8)
# }

# # legend(x=1, y=1.4, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "#494949", cex=2, pt.cex=4)
# dev.off()

# iso_overalldfkal <- iso_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal")
# # dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
# png(sprintf("./tx_noise_fig/diffiso_alltools_kal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
# ggplot(iso_overalldfkal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
#     geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
#     theme_light()+
#    scale_color_manual(values=colMap)+
#    facet_grid(rep~splits, labeller=seq_label)+
#    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=90, vjust = 1, hjust = 1), axis.title = element_text(size=20),
#         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
#         labs(color="Tools", shape="Background")
# dev.off()

# iso_overalldfsal <- iso_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# # dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
# png(sprintf("./tx_noise_fig/diffiso_alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
# ggplot(iso_overalldfsal, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
#     geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
#     theme_light()+
#    scale_color_manual(values=colMap)+
#    facet_grid(rep~splits, labeller=seq_label)+
#    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=90, vjust = 1, hjust = 1), axis.title = element_text(size=20),
#         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
#         labs(color="Tools", shape="Background")
# dev.off()


# png(sprintf("./tx_noise_fig/iso_alltools_%s_%s_%s_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
# ggplot(iso_overalldf, aes(x=recall, y=precision, color=tool, shape=as.character(noise)))+
#     geom_point(size=5, alpha=0.8)+geom_path(aes(x=recall, y=precision, group=tool))+
#     theme_light()+
#    scale_color_manual(values=colMap)+
#    facet_grid(con~splits, labeller=seq_label)+
#    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=90, vjust = 1, hjust = 1), axis.title = element_text(size=20),
#         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
#         labs(color="Tools", shape="Background")
# dev.off()


ev_overalldf <- get_overall(split="events")
ev_overalldf1 <- ev_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype)

#radar plot
library(fmsb)
library(gridExtra)

plot_list <- c()
plot_name <- c()
n<-1
for (rps in c("4 replicates", "8 replicates")){
    for (spl in sort(unique(ev_overalldf1$splits))){
        data <- ev_overalldf1 %>% filter(quant_tool=="salmon"&splits==spl&rep==rps) %>% select(noise, fdr, tool) %>% pivot_wider(names_from = tool, values_from = fdr, id_cols = noise) %>% as.data.frame()
        row.names(data) <- paste0("bg-",data$noise)
        data <- rbind(rep(1,ncol(data)) , rep(0,ncol(data)) , data)
        data <- data %>% select(-noise)
        data[is.na(data)] <- 0
        data <- lapply(data, function(x) if(is.factor(x) | is.character(x)) as.numeric(as.character(x)) else x) 

        plot_list[[n]] <- as.data.frame(data)
        plot_name[[n]] <- paste0(rps,"_",spl)
        n<-n+1
 }
}

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9), "#008000" , "black" )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4),"#008000" , "black" )

png(sprintf("./tx_noise_fig/ev_allradar_%s_%s_%s.png", seqtype, whichrep, thisrep),  width=3500, height=3000, res=300)
ev_overalldf %>%
    mutate(tool=fct_reorder(tool, fdr), fdr=as.numeric(fdr)) %>%
    filter(quant_tool=="rsem") %>%
ggplot(aes(x=tool, y=fdr, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=fdr, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(splits~con, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
        xlab("Tools")+ylab("FDR")+
        ggtitle("FDR in different event scenarios in paired-end data")
dev.off()
# png(sprintf("./tx_noise_fig/ev_allradar_%s_%s_%s.png", seqtype, whichrep, thisrep),  width=3500, height=2000, res=300)
# nrows<-2
# ncols<-3
# par(mfrow = c(nrows, ncols), mar = c(0.5,0.5,0.5,0.5), oma=c(4,4,4,4))  # 4 rows, 2 columns
#  # 4 rows, 2 columns

# for (i in seq_along(plot_list)) {
# radarchart(plot_list[[i]], axistype=1 , 
#         #custom polygon
#         pcol=colors_border ,  plwd=4 , plty=1,
#         #custom the grid
#         cglcol="grey", cglty=1, axislabcol="#5E5E5E", caxislabels=seq(0,0.6,0.15), cglwd=0.8,
#         #custom labels
#         vlcex=1.5, calcex=1
        
# )  
# }
# row_labels <- c("4 replicates", "8 replicates") # Adjust as needed
# for (i in 1:nrows) {
#   mtext(row_labels[i], side = 2, at = (nrows - i + 0.5)/nrows, srt=90, outer = TRUE, line = 1, cex=1.8)
# }

# # Add column labels
# col_labels <- c("S1", "S2", "S3")   # Adjust as needed
# for (i in 1:ncols) {
#   mtext(col_labels[i], side = 3, at = (i - 0.5)/ncols, las = 1, outer = TRUE, line = 1, cex=1.8)
# }

# # legend(x=1, y=1.4, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "#494949", cex=2, pt.cex=4)
# dev.off()

ev_overalldfsal <- ev_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./tx_noise_fig/ev_alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
ev_overalldfsal %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>% 
ggplot(aes(x=tool, y=fdr, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=fdr, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        labs(color="Tools", shape="Background")
dev.off()

png(sprintf("./tx_noise_fig/ev_alltools_%s_%s_%s_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
ev_overalldf %>%
    mutate(tool=fct_reorder(tool, desc(f1)), fdr=as.numeric(fdr)) %>%
ggplot(aes(x=tool, y=fdr, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=fdr, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(con~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
        ggtitle("FDR in different event scenarios in paired-end data")
dev.off()

png(sprintf("./tx_noise_fig/ev_alltools_%s_%s_%s_re_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
ev_overalldf %>%
    mutate(tool=fct_reorder(tool, desc(fdr))) %>%
ggplot(aes(x=tool, y=recall, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=recall, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(con~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
        ggtitle("Recall in different event scenarios in paired-end data")
dev.off()

png(sprintf("./tx_noise_fig/ev_alltools_%s_%s_%s_pr_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
ev_overalldf %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>%
ggplot(aes(x=tool, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(con~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
        ggtitle("Precision in different event scenarios in paired-end data")
dev.off()


gg_overalldf <- get_overall(split="gene_group")
gg_overalldf1 <- gg_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype)

#radar plot
library(fmsb)
library(gridExtra)

plot_list <- c()
plot_name <- c()
n<-1
for (rps in c("4 replicates", "8 replicates")){
    for (spl in sort(unique(gg_overalldf1$splits))){
        data <- gg_overalldf1 %>% filter(quant_tool=="salmon"&splits==spl&rep==rps) %>% select(noise, fdr, tool) %>% pivot_wider(names_from = tool, values_from = fdr, id_cols = noise) %>% as.data.frame()
        row.names(data) <- paste0("bg-",data$noise)
        data <- rbind(rep(1,ncol(data)) , rep(0,ncol(data)) , data)
        data <- data %>% select(-noise)
        data[is.na(data)] <- 0
        data <- lapply(data, function(x) if(is.factor(x) | is.character(x)) as.numeric(as.character(x)) else x)

        plot_list[[n]] <- as.data.frame(data)
        plot_name[[n]] <- paste0(rps,"_",spl)
        n<-n+1
 }
}

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9), "#008000" , "black" )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4),"#008000" , "black" )

png(sprintf("./tx_noise_fig/gg_allradar_%s_%s_%s.png", seqtype, whichrep, thisrep), width=3500, height=2000, res=300)
nrows<-2
ncols<-3
par(mfrow = c(nrows, ncols), mar = c(0.5,0.5,0.5,0.5), oma=c(4,4,4,4))  # 4 rows, 2 columns

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
col_labels <- sort(unique(gg_overalldf1$splits))   # Adjust as needed
for (i in 1:ncols) {
  mtext(col_labels[i], side = 3, at = (i - 0.5)/ncols, las = 1, outer = TRUE, line = 1, cex=1.8)
}

#legend(x=1, y=1.4, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "#494949", cex=2, pt.cex=4)
dev.off()


gg_overalldfkal <- gg_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="kal")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./tx_noise_fig/genegroup_alltools_kal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
gg_overalldfkal %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>% 
ggplot(aes(x=tool, y=f1, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=f1, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        labs(color="Tools", shape="Background")
dev.off()

gg_overalldfsal <- gg_overalldf %>% dplyr::filter(sim==whichrep&seq==seqtype&quant_tool=="salmon")
# dex <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/kal_res_gene_N_T.txt", sep="\t")
png(sprintf("./tx_noise_fig/genegroup_alltools_sal_%s_%s_%s.png", seqtype, whichrep, thisrep), width=900, height=600)
gg_overalldfsal %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>% 
ggplot(aes(x=tool, y=fdr, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=f1, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(rep~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        labs(color="Tools", shape="Background")
dev.off()


png(sprintf("./tx_noise_fig/gg_alltools_%s_%s_%s_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
gg_overalldf %>%
    mutate(tool=fct_reorder(tool, desc(f1)), fdr=as.numeric(fdr)) %>% 
ggplot(aes(x=tool, y=fdr, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=f1, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(con~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background.")+
        ggtitle("FDR scores in different gene group in paired-end data")
dev.off()

png(sprintf("./tx_noise_fig/gg_alltools_%s_%s_%s_re_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
gg_overalldf %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>% 
ggplot(aes(x=tool, y=recall, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=recall, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(con~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
        ggtitle("Recall in different gene group in paired-end data")
dev.off()

png(sprintf("./tx_noise_fig/gg_alltools_%s_%s_%s_pr_supp.png", seqtype, whichrep, thisrep), width=800, height=1200)
gg_overalldf %>%
    mutate(tool=fct_reorder(tool, desc(f1))) %>% 
ggplot(aes(x=tool, y=precision, color=tool, shape=as.character(noise)))+
    geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=precision, group=tool))+
    theme_light()+
   scale_color_manual(values=colMap)+
   facet_grid(con~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
        labs(color="Tools", shape="Background")+
        ggtitle("Precision in different gene group in paired-end data")
dev.off()
# tru <- tru %>% dplyr::filter(fc_group==4)
# p_g <- tru$feature_id[tru$status==1]
# dex1 <- dex %>% dplyr::mutate(dexseq=replace_na(dexseq, 1))
# dex1 <- dex1 %>% group_by(feature_id) %>% summarise(across(colnames(select(dex, -feature_id)), .fns = min), .groups = "keep")

# d_g <- dex1$feature_id[dex1$dexseq<0.05]

# sal <- read.csv("/nfs/scratch/chit/simulated_real/single_50_8_r1/results/rsem_count.csv")
# sum(is.na(dex$feature_id))
View(as.data.frame(fc_overalldf)) 
tru <- read.csv("/nfs/scratch/chit/simulated_real/pair_50_8_r1_0.5/results/truthtable_gene.csv", sep="\t")

diffiso <- read.csv("/nfs/scratch/chit/simulated_real/pair_50_8_r1_0.5/results/isoinfo.csv", sep="\t")
isorat <- diffiso %>% group_by(gene_id) %>% summarise(diffiso=max(diffiso))
tru$diffiso <- isorat$diffiso
tru <- tru %>% dplyr::mutate(diffiso=ifelse(diffiso>0.8, ">0.8", ifelse(diffiso>=0.6, "0.6-0.8", ifelse(diffiso>=0.5, "0.5-0.6", ifelse(diffiso>=0.4, "0.4-0.5", ifelse(diffiso>=0.3, "0.3-0.4", ifelse(diffiso>=0.2, "0.2-0.3", ifelse(diffiso>=0.1, "0.1-0.2", "0-0.1")))))))) 


png("./tx_noise_fig/tru_stat.png")
ggplot(tru, aes(x=diffiso))+
    geom_bar()+
    facet_grid(gene_group~fc_group)
dev.off()
# sum(d_g %in% p_g)/length(p_g)
sal <- read.csv("/nfs/scratch/chit/simulated_real/pair_50_8_r1_0.5/results/salmon_res_tx_N_T.txt", sep="\t")

sal <- sal[grepl("ENSG", sal$feature_id),]


fc_t <- tru %>% dplyr::filter(events=="none")
fcsplitted <- sal %>% filter(feature_id %in% fc_t$feature_id)
fcval <- lapply(fcsplitted$iso_ktsp, function(x){ifelse(is.na(x), 0, x)}) %>% unlist
sum(fcval > 0.5)

gg_t <- tru %>% dplyr::filter(gene_group == ">9") #%>% filter(status==1)
ggsplitted <- sal %>% filter(feature_id %in% gg_t$feature_id)
ggval <- lapply(ggsplitted$iso_ktsp, function(x){ifelse(is.na(x), 0, x)}) %>% unlist
sum(ggval > 0.5)

tru %>% filter(fc_group==1)
fil1 <- sal[sal$feature_id %in% tru$feature_id,]

detectedpos <- fil1$feature_id[fil1$drimseq <0.05]
boomet <- unique(detectedpos)
tru <- tru %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
 
split_tru <- tru %>% filter(status==1)
split_tru %>% nrow
detectedpos[detectedpos %in% split_tru$feature_id] 

tru %>% head

fc_overalldf %>% filter(tool=="nbsplice" & quant_tool=="salmon")


nrow(tru %>% dplyr::filter(fc_group=="2"&status==1))+
    nrow(tru %>% dplyr::filter(fc_group=="3"&status==1))+
    nrow(tru %>% dplyr::filter(fc_group=="4"&status==1))+
    length(tru %>% dplyr::filter(fc_group=="5"&status==1))


3975+4118+3877+4042

nrow(tru %>% dplyr::filter(gene_group=="2-4"&status==1))+
    nrow(tru %>% dplyr::filter(gene_group=="5-9"&status==1))+
    nrow(tru %>% dplyr::filter(gene_group==">9"&status==1))
