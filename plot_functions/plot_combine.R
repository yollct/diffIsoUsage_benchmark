library(tidyverse)
library(ggplot2)
library(UpSetR)
library(ggsci)
library(RColorBrewer)
source("/nfs/home/students/chit/is_benchmark/plot_functions/functions.R")

path <- "/nfs/home/students/chit/is_benchmark"
outdir <- "/nfs/scratch/chit/simulated_real/%s_50_%s_%s"

get_overall <- function(split){
    overalldf <- do.call(rbind, lapply(c("r2"), function(R){
            do.call(rbind, lapply(c("single", "pair"), function(seq){
            do.call(rbind, lapply(c("8"), function(rep){
                truthfiles_g <- read.csv(paste0(sprintf(outdir, seq, rep, R), "/results/truthtable_gene.csv"), sep="\t") ### change to result
                truthfiles_tx <- read.csv(paste0(sprintf(outdir, seq, rep, R), "/results/truthtable_tx.csv"), sep="\t")

                # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                # truthfiles_tx <- truthfiles_tx %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                
                do.call(rbind, lapply(c("salmon", "kal", "rsem"), function(quant){
                    df_g <- read.csv(paste0(sprintf(outdir, seq,rep, R), "/results/", sprintf("%s_res_gene_N_T.txt", quant)), sep="\t")
                    thisthres <- thresholds[colnames(df_g)[2:ncol(df_g)],]$thres

                    outputpr <- cal_pre_re(df_g, truthfiles_g, split=split)
                    this_pr <- pivot_output(outputpr, split=split)
                    this_pr$quant_tool <- quant
                    this_pr$rep <- rep
                    this_pr$seq <- seq
                    this_pr$sim <- R
                    rm(df_g)
                    this_pr
                }))
            }))
        }))
    }))
    overalldf
}
seq_label <- as_labeller(c('pair'="Paired-end", 'single'='Single-end', '2'='2', '3'='3','4'='4','5'='5', '2-4'='2-4','5-9'='5-9','>9'='>9','DTE'='DTE','DTU'='DTU','DTEDTU'='DTEDTU'))

overalldf <- get_overall(split=NULL)
overalldf1 <- overalldf %>% dplyr::filter(sim=="r2")
png("alltools_pr_8.png", width=800, height=400)
ggplot(overalldf1, aes(x=precision, y=recall, color=tool, shape=quant_tool)) +
  geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_nejm()+
  facet_grid(.~seq, labeller=seq_label)+
    theme(axis.text=element_text(size=20), axis.text.x = element_text(vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Precision")+
    ylab("Recall")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()

### gene group
fc_overalldf <- get_overall(split="fc_group")
fc_overalldf1 <- fc_overalldf %>% dplyr::filter(sim=="r2")
png("fc_alltools_pr.png", width=900, height=600)
ggplot(fc_overalldf1, aes(x=precision, y=recall, color=tool, shape=quant_tool))+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_nejm()+
   facet_grid(seq~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Precision")+
    ylab("Recall")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()

iso_overalldf <- get_overall(split="gene_group")
iso_overalldf1 <- iso_overalldf %>% dplyr::filter(sim=="r2")
png("iso_alltools_pr.png", width=900, height=600)
ggplot(iso_overalldf1, aes(x=precision, y=recall, color=tool, shape=quant_tool))+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_nejm()+
   facet_grid(seq~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Precision")+
    ylab("Recall")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()

ev_overalldf <- get_overall(split="events")
ev_overalldf1 <- ev_overalldf %>% dplyr::filter(sim=="r1")
png("ev_alltools_pr_8.png", width=900, height=600)
ggplot(ev_overalldf1, aes(x=precision, y=recall, color=tool, shape=quant_tool))+
    geom_point(size=7, alpha=0.8)+ theme_light()+
   scale_color_nejm()+
   facet_grid(seq~splits, labeller=seq_label)+
   theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
    xlab("Precision")+
    ylab("Recall")+
    scale_shape_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    labs(color="Tools", shape="Quantification tools")
dev.off()

mean_quant <- overalldf %>% group_by(tool, quant_tool, rep, seq) %>% summarise(nmin_p = min(precision), nmax_p = max(precision), nmin_r = min(recall), nmax_r=max(recall), precision=mean(precision), recall=mean(recall))
png("alltools_p.png", width=500, height=500)
ggplot(mean_quant, aes(x=quant_tool, y=precision, color=tool, shape=tool, ymin=nmin_r, ymax=nmax_r)) +

  geom_point(size=5, stroke=1.5, alpha=0.7)+ theme_light()+
  scale_shape_manual(values=1:nlevels(overalldf$tool))+
  facet_grid(rep~seq)+
    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))
dev.off()
### only tools with diff quant tools
#overall_quant <- overalldf %>% dplyr::filter(tool!="junctionseq", tool!="seqGSEA" , tool!="cuffdiff")


### for quant corr
overallquant <- do.call(rbind, lapply(c("r1", "r2", "r3"), function(R){
    do.call(rbind, lapply(c("single", "pair"), function(seq){
        do.call(rbind, lapply(c("4", "8"), function(rep){
            groundtruthfiles <- list.files(paste0(sprintf(outdir, seq, rep, R), "/fastq_sim"), pattern="isoforms.results")### change to result
            groundtruth <- do.call("cbind", lapply(groundtruthfiles, function(x){
                        tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R),"/fastq_sim/", x), sep="\t")
                        tmp$count
            }))
            colnames(groundtruth) <- lapply(groundtruthfiles, function(x){gsub(".sim.isoforms.results", "", x)})
            tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R),"/fastq_sim/", groundtruthfiles[1]), sep="\t")
            groundtruth <- as.data.frame(groundtruth)
            groundtruth$tx <- tmp$transcript_id
            
            groundtruth_long <- pivot_longer(groundtruth, cols=starts_with("SRR"), names_to="sample", values_to = "count")

            do.call(rbind, lapply(c("salmon_out", "kallisto_out", "rsem_out"), function(quant){
                print(quant)
                quant_pattern <- ifelse(quant=="salmon_out", "quant.sf", ifelse(quant=="kallisto_out", "abundance.tsv", ".isoforms.results"))
                txcol <- ifelse(quant=="salmon_out", "Name", ifelse(quant=="kallisto_out", "target_id", "transcript_id"))

                res_files <- list.files(paste0(sprintf(outdir, seq,rep,R),"/",quant))
                df_g <- do.call("cbind", lapply(res_files, function(x){
                    quantfiles <- ifelse(quant=="rsem_out", sprintf("%s%s",x,quant_pattern), sprintf("%s", quant_pattern))
                    quantcol <- ifelse(quant=="rsem_out", "expected_count", ifelse(quant=="salmon_out", "NumReads", "est_counts"))

                    tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R), "/", quant, "/", x, "/" ,quantfiles), sep="\t")
                    tmp[,quantcol]
                    
                }))
                df_g <- as.data.frame(df_g)
                quanttmp <- ifelse(quant=="rsem_out", sprintf("%s/%s%s",res_files[1],res_files[1],quant_pattern), sprintf("%s/%s", res_files[1], quant_pattern))

                tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R), "/", quant, "/" ,quanttmp), sep="\t")
                colnames(df_g) <- res_files
                df_g$tx <- lapply(tmp[,txcol], function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist

                quant_long <- pivot_longer(df_g, cols=starts_with("SRR"), names_to="sample", values_to="count")

                df_gt <- full_join(groundtruth_long, quant_long, by=c("tx","sample"))
                df_gt <- df_gt %>% na.omit()

                corr_res <- data.frame('spearman'=cor(df_gt$count.x, df_gt$count.y, method = "spearman"),
                                        'pearson'=cor(df_gt$count.x, df_gt$count.y, method = "pearson"))

                corr_res$quant_tool <- gsub("_out", "", quant)
                corr_res$seq <- seq
                corr_res$rep <- rep
                corr_res$sim <- R

                rm(df_g)
                return(corr_res)
            }))
        }))
    }))
}))


png("overquant_pr.png", width=500, height=500)
ggplot(overalldf %>% dplyr::filter(tool!="iso_ktsp") %>% dplyr::filter(tool!="junctionseq") %>% dplyr::filter(tool!="cuffdiff") %>% dplyr::filter(tool!="seqGSEA"), aes(x=quant_tool, y=precision, fill=quant_tool, shape=quant_tool))+
    geom_violin()+
    scale_fill_npg(name="Quantification", palette = c("nrc"))+
    scale_shape_manual(values=1:nlevels(overalldf$tool))+
    scale_x_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    ylim(c(0.4,1))+
    xlab("Quantification tool")+
    ylab("Precision")+
    facet_grid(.~seq, labeller=seq_label)+ theme_light()+
    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        guides(fill="none")

dev.off()



png("overquant_re.png", width=500, height=500)
ggplot(overalldf %>% dplyr::filter(tool!="iso_ktsp") %>% dplyr::filter(tool!="junctionseq") %>% dplyr::filter(tool!="cuffdiff") %>% dplyr::filter(tool!="seqGSEA"), aes(x=quant_tool, y=as.numeric(recall), fill=quant_tool, shape=quant_tool))+
    geom_violin()+
    scale_shape_manual(values=1:nlevels(overalldf$tool))+
    scale_x_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))+
    xlab("Quantification tool")+
    ylab("Recall")+
    ylim(0,0.6)+
    facet_grid(.~seq, labeller=seq_label)+ theme_light()+
    theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=45, vjust = 1, hjust = 1), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        scale_fill_npg(name="Quantification", palette = c("nrc"))+
        guides(fill = "none")
dev.off()

png("overquant.png", width=500, height=500)
ggplot(overallquant, aes(y=spearman, x=seq, color=quant_tool))+
    geom_boxplot()+
    theme_light()+theme(axis.text=element_text(size=20), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        scale_color_npg(name="Quantification", palette = c("nrc"))+
        xlab("Sequencing type")+scale_x_discrete(label = c("Paired-end", "Single-end"))
dev.off()        

long_quant <- pivot_longer(overallquant, cols=c("spearman", "pearson"), names_to = "corr", values_to="value", )
png("overquant_sp.png", width=1000, height=500)
ggplot(long_quant, aes(y=value, x=seq, color=quant_tool))+
    geom_boxplot()+
    facet_grid(.~corr)+scale_color_npg(name="Quantification", palette = c("nrc"))+
    theme_light()+theme(axis.text=element_text(size=20), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        xlab("Sequencing type")+scale_x_discrete(label = c("Paired-end", "Single-end"))+ylab("Correlation")
dev.off()

### count how many expressed isoforms are in each condition
count_express <- do.call(rbind, lapply(c("r1", "r2", "r3"), function(R){
    do.call(rbind, lapply(c("single", "pair"), function(seq){
        do.call(rbind, lapply(c("4", "8"), function(rep){
            groundtruthfiles <- list.files(paste0(sprintf(outdir, seq, rep, R), "/fastq_sim"), pattern="isoforms.results")### change to result
            groundtruth <- do.call("cbind", lapply(groundtruthfiles, function(x){
                        tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R),"/fastq_sim/", x), sep="\t")
                        tmp$count
            }))
            colnames(groundtruth) <- lapply(groundtruthfiles, function(x){gsub(".sim.isoforms.results", "", x)})
            tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R),"/fastq_sim/", groundtruthfiles[1]), sep="\t")
            groundtruth <- as.data.frame(groundtruth)
            groundtruth$tx <- tmp$transcript_id
            
            groundtruth_long <- pivot_longer(groundtruth, cols=starts_with("SRR"), names_to="sample", values_to = "count")

            do.call(rbind, lapply(c("salmon_out", "kallisto_out", "rsem_out"), function(quant){
                print(quant)
                quant_pattern <- ifelse(quant=="salmon_out", "quant.sf", ifelse(quant=="kallisto_out", "abundance.tsv", ".isoforms.results"))
                txcol <- ifelse(quant=="salmon_out", "Name", ifelse(quant=="kallisto_out", "target_id", "transcript_id"))
                quantcol <- ifelse(quant=="rsem_out", "expected_count", ifelse(quant=="salmon_out", "NumReads", "est_counts"))
                res_files <- list.files(paste0(sprintf(outdir, seq,rep,R),"/",quant))

                df_g <- do.call("cbind", lapply(res_files, function(x){
                    quantfiles <- ifelse(quant=="rsem_out", sprintf("%s%s",x,quant_pattern), sprintf("%s", quant_pattern))
                    
                    tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R), "/", quant, "/", x, "/" ,quantfiles), sep="\t")
                    tmp[,quantcol]
                    
                }))
                df_g <- as.data.frame(df_g)
                quanttmp <- ifelse(quant=="rsem_out", sprintf("%s/%s%s",res_files[1],res_files[1],quant_pattern), sprintf("%s/%s", res_files[1], quant_pattern))

                tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R), "/", quant, "/" ,quanttmp), sep="\t")
                colnames(df_g) <- res_files
                df_g$tx <- lapply(tmp[,txcol], function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
                
                n_expresssed<-sum(tmp[,quantcol]>1)
                
                corr_res = data.frame(n_expresssed=n_expresssed, seq=seq, rep=rep, tool=quant)

                rm(df_g)
                return(corr_res)
            }))
        }))
    }))
}))

png("n_transcript.png", width=500, height=500)
ggplot(count_express, aes(x=seq, y=n_expresssed, fill=tool)) +
    geom_violin()+
    theme_light()+theme(axis.text=element_text(size=20), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        scale_fill_npg(name="Quantification", palette = c("nrc"))+
        xlab("Sequencing type")+scale_x_discrete(label = c("Paired-end", "Single-end"))+
        ylab("Number of expressed isoforms (count > 10)")
dev.off()


#### RMSE
quantrmse <- do.call(rbind, lapply(c("r1", "r2", "r3"), function(R){
    do.call(rbind, lapply(c("single", "pair"), function(seq){
        do.call(rbind, lapply(c("4", "8"), function(rep){
            groundtruthfiles <- list.files(paste0(sprintf(outdir, seq, rep, R), "/fastq_sim"), pattern="isoforms.results")### change to result
            groundtruth <- do.call("cbind", lapply(groundtruthfiles, function(x){
                        tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R),"/fastq_sim/", x), sep="\t")
                        tmp$count
            }))
            colnames(groundtruth) <- lapply(groundtruthfiles, function(x){gsub(".sim.isoforms.results", "", x)})
            tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R),"/fastq_sim/", groundtruthfiles[1]), sep="\t")
            groundtruth <- as.data.frame(groundtruth)
            groundtruth$tx <- tmp$transcript_id
            
            groundtruth_long <- pivot_longer(groundtruth, cols=starts_with("SRR"), names_to="sample", values_to = "count")

            do.call(rbind, lapply(c("salmon_out", "kallisto_out", "rsem_out"), function(quant){
                print(quant)
                quant_pattern <- ifelse(quant=="salmon_out", "quant.sf", ifelse(quant=="kallisto_out", "abundance.tsv", ".isoforms.results"))
                txcol <- ifelse(quant=="salmon_out", "Name", ifelse(quant=="kallisto_out", "target_id", "transcript_id"))

                res_files <- list.files(paste0(sprintf(outdir, seq,rep,R),"/",quant))
                df_g <- do.call("cbind", lapply(res_files, function(x){
                    quantfiles <- ifelse(quant=="rsem_out", sprintf("%s%s",x,quant_pattern), sprintf("%s", quant_pattern))
                    quantcol <- ifelse(quant=="rsem_out", "expected_count", ifelse(quant=="salmon_out", "NumReads", "est_counts"))

                    tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R), "/", quant, "/", x, "/" ,quantfiles), sep="\t")
                    tmp[,quantcol]
                    
                }))
                df_g <- as.data.frame(df_g)
                quanttmp <- ifelse(quant=="rsem_out", sprintf("%s/%s%s",res_files[1],res_files[1],quant_pattern), sprintf("%s/%s", res_files[1], quant_pattern))

                tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R), "/", quant, "/" ,quanttmp), sep="\t")
                colnames(df_g) <- res_files
                df_g$tx <- lapply(tmp[,txcol], function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist

                quant_long <- pivot_longer(df_g, cols=starts_with("SRR"), names_to="sample", values_to="count")

                df_gt <- full_join(groundtruth_long, quant_long, by=c("tx","sample"))
                df_gt <- df_gt %>% na.omit()

                corr_res <- data.frame('rmse'= ((df_gt$count.y-df_gt$count.x)**2))

                corr_res$quant_tool <- gsub("_out", "", quant)
                corr_res$seq <- seq
                corr_res$rep <- rep
                corr_res$sim <- R

                rm(df_g)
                return(corr_res)
            }))
            }))
        }))
    }))

overallrmse <- quantrmse %>% dplyr::filter(rmse!=0) %>% group_by(quant_tool, seq, rep, sim) %>% summarise(rmse=sqrt(mean(rmse)))

png("rmse_quant.png", width=900, height=500)
ggplot(overallrmse, aes(x=quant_tool, y=rmse, fill=quant_tool))+
    geom_boxplot()+
    facet_wrap(.~seq, labeller=seq_label, scales="free_y")+
    theme_light()+theme(axis.text=element_text(size=20), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        scale_fill_npg(name="Quantification", palette = c("nrc"))+
        xlab("Quantification")+
        ylab("RMSE")
dev.off()

