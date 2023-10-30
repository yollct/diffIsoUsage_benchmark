library(tidyverse)
library(ggplot2)
library(UpSetR)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(GGally)
source("/nfs/home/students/chit/is_benchmark/plot_functions/functions.R")

path <- "/nfs/home/students/chit/is_benchmark"
outdir <- "/nfs/scratch/chit/simulated_real/%s_50_%s_%s"

thresholds <- data.frame(tools=c("iso_ktsp","drimseq","dexseq","dturtle","seqGSEA","cuffdiff","junctionseq","saturn","drimseq_stageR", "dexseq_stageR","saturn_stageR", "DSGseq", "nbsplice", "LimmaDS", "edgeR"),
                        thres=c(0.8,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05, 2, 0.05, 0.05, 0.05))
row.names(thresholds) <- thresholds$tools

noise <- ""
thisrep <- c('4')
seqtype <- 'single'
whichrep <- 'r1'
get_overall1 <- function(split){
    overalldf <- do.call(rbind, lapply(whichrep, function(R){
            do.call(rbind, lapply(seqtype, function(seq){
            do.call(rbind, lapply(thisrep, function(rep){
                #do.call(rbind, lapply(noises, function(noise){
                    print(paste0(sprintf(outdir, seq, rep, R)))
                    truthfiles_g <- read.csv(paste0(sprintf(outdir, seq, rep, R), "/results/truthtable_gene.csv"), sep="\t") ### change to result
                    truthfiles_tx <- read.csv(paste0(sprintf(outdir, seq, rep, R), "/results/truthtable_tx.csv"), sep="\t")
                    
                    
                    # if (is.null(split)){
                    #     print("not consider DTE")
                    #     truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                    # } else {
                    #     if (split!="events"){
                    #     print("not consider DTE")
                    #     truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                    #     }
                    # }

                    diffiso <- read.csv(paste0(sprintf(outdir, seq, rep, R), "/results/isoinfo.csv"), sep="\t")
                    isorat <- diffiso %>% group_by(gene_id) %>% summarise(diffiso=max(diffiso))
                    truthfiles_g$diffiso <- isorat$diffiso
                    truthfiles_g <- truthfiles_g %>% dplyr::mutate(diffiso=ifelse(diffiso>0.8, ">0.8", ifelse(diffiso>=0.6, "0.6-0.8", ifelse(diffiso>=0.5, "0.5-0.6", ifelse(diffiso>=0.4, "0.4-0.5", ifelse(diffiso>=0.3, "0.3-0.4", ifelse(diffiso>=0.2, "0.2-0.3", ifelse(diffiso>=0.1, "0.1-0.2", "0-0.1")))))))) 
                    print(truthfiles_g %>% nrow)
                    # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    # truthfiles_tx <- truthfiles_tx %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    
                    do.call(rbind, lapply(c("salmon", "kal", "rsem"), function(quant){
                        df_g <- read.csv(paste0(sprintf(outdir, seq,rep, R), "/results/", sprintf("%s_res_gene_N_T.txt", quant)), sep="\t")
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
 
    overalldf$f1 <- lapply(overalldf$f1, function(x){as.numeric(x)}) %>% unlist
    overalldf
}

seq_label <- as_labeller(c('pair'="Paired-end", 'single'='Single-end', '2'='2', '3'='3','4'='4','5'='5', '2-4'='2-4','5-9'='5-9','>9'='>9','DTE'='DTE','DTU'='DTU','DTEDTU'='DTEDTU','0-0.3'='0-0.3','>0.8'='>0.8', '0.4-0.7'='0.4-0.7', '8 replicates'='8 replicates', '4 replicates'='4 replicates', "_0.1"="Noise 0.1", "_0.5"="Noise 0.5"))


overalldf <- get_overall1(split=NULL)

### for quant corr
overallquant <- do.call(rbind, lapply(c("", "_0.1", "_0.5"), function(noise){
    do.call(rbind, lapply(c("r1"), function(R){
        do.call(rbind, lapply(c("single", "pair"), function(seq){
            do.call(rbind, lapply(c("4", "8"), function(rep){
                groundtruthfiles <- list.files(paste0(sprintf(outdir, seq, rep, R, noise), "/fastq_sim"), pattern="isoforms.results")### change to result
                groundtruth <- do.call("cbind", lapply(groundtruthfiles, function(x){
                            tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R, noise),"/fastq_sim/", x), sep="\t")
                            tmp$count
                }))
                colnames(groundtruth) <- lapply(groundtruthfiles, function(x){gsub(".sim.isoforms.results", "", x)})
                tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R, noise),"/fastq_sim/", groundtruthfiles[1]), sep="\t")
                groundtruth <- as.data.frame(groundtruth)
                groundtruth$tx <- tmp$transcript_id
                
                groundtruth_long <- pivot_longer(groundtruth, cols=starts_with("SRR"), names_to="sample", values_to = "count")

                do.call(rbind, lapply(c("salmon_out", "kallisto_out", "rsem_out"), function(quant){
                    print(quant)
                    quant_pattern <- ifelse(quant=="salmon_out", "quant.sf", ifelse(quant=="kallisto_out", "abundance.tsv", ".isoforms.results"))
                    txcol <- ifelse(quant=="salmon_out", "Name", ifelse(quant=="kallisto_out", "target_id", "transcript_id"))

                    res_files <- list.files(paste0(sprintf(outdir, seq,rep,R, noise),"/",quant))
                    df_g <- do.call("cbind", lapply(res_files, function(x){
                        quantfiles <- ifelse(quant=="rsem_out", sprintf("%s%s",x,quant_pattern), sprintf("%s", quant_pattern))
                        quantcol <- ifelse(quant=="rsem_out", "expected_count", ifelse(quant=="salmon_out", "NumReads", "est_counts"))

                        tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R, noise), "/", quant, "/", x, "/" ,quantfiles), sep="\t")
                        tmp[,quantcol]
                        
                    }))
                    df_g <- as.data.frame(df_g)
                    quanttmp <- ifelse(quant=="rsem_out", sprintf("%s/%s%s",res_files[1],res_files[1],quant_pattern), sprintf("%s/%s", res_files[1], quant_pattern))

                    tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R, noise), "/", quant, "/" ,quanttmp), sep="\t")
                    colnames(df_g) <- res_files
                    df_g$tx <- lapply(tmp[,txcol], function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist

                    quant_long <- pivot_longer(df_g, cols=starts_with("SRR"), names_to="sample", values_to="count")

                    df_gt <- full_join(groundtruth_long, quant_long, by=c("tx","sample"))
                    df_gt <- df_gt %>% na.omit()
                    print(df_gt %>% head)

                    corr_res <- data.frame('spearman'=cor(df_gt$count.x, df_gt$count.y, method = "spearman"),
                                            'pearson'=cor(df_gt$count.x, df_gt$count.y, method = "pearson"),
                                            'quant_tool'=gsub("_out", "", quant),
                                            'seq'=seq,
                                            'rep'=rep,
                                            'sim'=R,
                                            'tx'=df_gt$tx,
                                            'noise' = noise)

                    
                    rm(df_g)
                    return(corr_res)
                }))
            }))
        }))
    }))
}))




png("overquant_pear.png", width=800, height=500)
ggplot(overallquant, aes(y=pearson, x=quant_tool))+
    geom_boxplot(fill="#5cbebb")+
    facet_grid(rep~seq,labeller=seq_label, scales="free")+
    ylim(c(0,1))+
    theme_light()+theme(axis.text=element_text(size=20), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        scale_color_npg(name="Quantification", palette = c("nrc"))+
        xlab("Quantification")+scale_x_discrete(labels=c("kal" = "Kallisto", "rsem" = "RSEM",
                              "salmon" = "Salmon"))
dev.off()        

long_quant <- pivot_longer(overallquant, cols=c("spearman", "pearson"), names_to = "corr", values_to="value", )
png("overquant_sppr.png", width=1000, height=500)
ggplot(long_quant, aes(y=value, x=seq, color=quant_tool))+
    geom_boxplot()+
    facet_grid(rep~corr)+scale_color_d3(name="Quantification", palette = c("category20"))+
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
        scale_fill_d3(name="Quantification", palette = c("category20"))+
        xlab("Sequencing type")+scale_x_discrete(label = c("Paired-end", "Single-end"))+
        ylab("Number of expressed isoforms (count > 10)")
dev.off()


#### RMSE
quantrmse <- do.call(rbind, lapply(c("", "_0.1", "0.5"), function(noise){
        do.call(rbind, lapply(c("r1"), function(R){
            do.call(rbind, lapply(c("single", "pair"), function(seq){
                do.call(rbind, lapply(c("4", "8"), function(rep){
                    groundtruthfiles <- list.files(paste0(sprintf(outdir, seq, rep, R, noise), "/fastq_sim"), pattern="isoforms.results")### change to result
                    groundtruth <- do.call("cbind", lapply(groundtruthfiles, function(x){
                                tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R, noise),"/fastq_sim/", x), sep="\t")
                                tmp$count
                    }))
                    colnames(groundtruth) <- lapply(groundtruthfiles, function(x){gsub(".sim.isoforms.results", "", x)})
                    tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R, noise),"/fastq_sim/", groundtruthfiles[1]), sep="\t")
                    groundtruth <- as.data.frame(groundtruth)
                    groundtruth$tx <- tmp$transcript_id
                    
                    groundtruth_long <- pivot_longer(groundtruth, cols=starts_with("SRR"), names_to="sample", values_to = "count")

                    do.call(rbind, lapply(c("salmon_out", "kallisto_out", "rsem_out"), function(quant){
                        print(quant)
                        quant_pattern <- ifelse(quant=="salmon_out", "quant.sf", ifelse(quant=="kallisto_out", "abundance.tsv", ".isoforms.results"))
                        txcol <- ifelse(quant=="salmon_out", "Name", ifelse(quant=="kallisto_out", "target_id", "transcript_id"))

                        res_files <- list.files(paste0(sprintf(outdir, seq,rep,R, noise),"/",quant))
                        df_g <- do.call("cbind", lapply(res_files, function(x){
                            quantfiles <- ifelse(quant=="rsem_out", sprintf("%s%s",x,quant_pattern), sprintf("%s", quant_pattern))
                            quantcol <- ifelse(quant=="rsem_out", "expected_count", ifelse(quant=="salmon_out", "NumReads", "est_counts"))

                            tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R, noise), "/", quant, "/", x, "/" ,quantfiles), sep="\t")
                        
                            thistmp <- tmp[,quantcol]
                            thistmp
                        }))
                        df_g <- as.data.frame(df_g)
                        
                        quanttmp <- ifelse(quant=="rsem_out", sprintf("%s/%s%s",res_files[1],res_files[1],quant_pattern), sprintf("%s/%s", res_files[1], quant_pattern))

                        tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R, noise), "/", quant, "/" ,quanttmp), sep="\t")
                        
                        colnames(df_g) <- res_files
                        df_g$tx <- lapply(tmp[,txcol], function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
                        
                        quant_long <- pivot_longer(df_g, cols=starts_with("SRR"), names_to="sample", values_to="count")
                        
                        df_gt <- full_join(groundtruth_long, quant_long, by=c("tx","sample"))
                        df_gt <- df_gt %>% na.omit()
            
                        corr_res <- data.frame('rmse'= sqrt(((df_gt$count.y-df_gt$count.x)**2)/nrow(df_gt)), 
                                                'quant_tool'=gsub("_out", "", quant),
                                                'seq'=seq,
                                                'rep'=rep,
                                                'sim'=R, 
                                                'tx'=df_gt$tx,
                                                'noise'=noise)
                        
                        rm(df_g)
                        return(corr_res)
                    }))
                }))
            }))
        }))
}))

overallrmse <- quantrmse %>% dplyr::filter(rmse!=0) %>% group_by(quant_tool, seq, rep, sim, noise) %>% summarise(rmse=sqrt(mean(rmse)))

png("rmse_quant.png", width=900, height=500)
ggplot(overallrmse, aes(x=quant_tool, y=rmse))+
    geom_boxplot(fill="#5cbebb")+
    facet_wrap(rep~seq, labeller=seq_label, scales="free_y")+
    theme_light()+theme(axis.text=element_text(size=20), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        scale_fill_d3(name="Quantification", palette = c("category20"))+
        xlab("Quantification")+
        ylab("RMSE")
dev.off()


####### get F1 score ########
get_f1_scores <- function(split){
    overalldf <- do.call(rbind, lapply(c("r1"), function(R){
            do.call(rbind, lapply(c("single", "pair"), function(seq){
            do.call(rbind, lapply(c("8"), function(rep){
                truthfiles_g <- read.csv(paste0(sprintf(outdir, seq, rep, R), "/results/truthtable_gene.csv"), sep="\t") ### change to result
                truthfiles_tx <- read.csv(paste0(sprintf(outdir, seq, rep, R), "/results/truthtable_tx.csv"), sep="\t")
                
                if (!is.null(split)&&split!="events" ){
                    truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                }
                diffiso <- read.csv(paste0(sprintf(outdir, 'single', '4', 'r1'), "/results/isoinfo.csv"), sep="\t")
                isorat <- diffiso %>% group_by(gene_id) %>% summarise(diffiso=max(diffiso))
                truthfiles_g$diffiso <- isorat$diffiso
                truthfiles_g <- truthfiles_g %>% dplyr::mutate(diffiso=ifelse(diffiso>0.8, ">0.8", ifelse(diffiso>0.4, "0.4-0.7", "0-0.3")))
                # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                # truthfiles_tx <- truthfiles_tx %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                
                do.call(rbind, lapply(c("salmon", "kal", "rsem"), function(quant){
                    df_g <- read.csv(paste0(sprintf(outdir, seq,rep, R), "/results/", sprintf("%s_res_gene_N_T.txt", quant)), sep="\t")
                    thisthres <- thresholds[colnames(df_g)[2:ncol(df_g)],]$thres

                    outputpr <- cal_f1(df_g, truthfiles_g, split=split)
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

f1df <- get_f1_scores(split=NULL)

rmsemin <- min(quantrmse$rmse)
rmsemax <- max(quantrmse$rmse)
quantrmse <- quantrmse %>% dplyr::mutate(rmse_norm=(rmse-rmsemin)/(rmsemax-rmsemin))
quantrmse %>% head

overallquant$forjoin <- paste0(overallquant$tx,overallquant$quant_tool,overallquant$seq, overallquant$rep, overallquant$sim)
quantrmse$forjoin <- paste0(quantrmse$tx, quantrmse$quant_tool, quantrmse$seq, quantrmse$rep, quantrmse$sim)
allquant <- left_join(overallquant, quantrmse, by=c("forjoin"))

png("/nfs/home/students/chit/is_benchmark/noise_fig/tx_compare.png")
ggplot(allquant, aes(x=spearman, y=rmse_norm))+
    geom_point()+
    facet_grid(quant_tool.x~seq.y)
dev.off()

png("/nfs/home/students/chit/is_benchmark/noise_fig/tx_compare_pearson.png")
ggplot(allquant, aes(x=pearson, y=rmse_norm))+
    geom_point()+
    facet_grid(quant_tool.x~seq.y)
dev.off()


#### RMSE
quantaberr <- do.call(rbind, lapply(c("", "_0.1", "0.5"), function(noise){
        do.call(rbind, lapply(c("r1"), function(R){
            do.call(rbind, lapply(c("single", "pair"), function(seq){
                do.call(rbind, lapply(c("4", "8"), function(rep){
                    groundtruthfiles <- list.files(paste0(sprintf(outdir, seq, rep, R, noise), "/fastq_sim"), pattern="isoforms.results")### change to result
                    groundtruth <- do.call("cbind", lapply(groundtruthfiles, function(x){
                                tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R, noise),"/fastq_sim/", x), sep="\t")
                                tmp$count
                    }))
                    colnames(groundtruth) <- lapply(groundtruthfiles, function(x){gsub(".sim.isoforms.results", "", x)})
                    tmp<-read.csv(paste0(sprintf(outdir, seq, rep, R, noise),"/fastq_sim/", groundtruthfiles[1]), sep="\t")
                    groundtruth <- as.data.frame(groundtruth)
                    groundtruth$tx <- tmp$transcript_id
                    
                    groundtruth_long <- pivot_longer(groundtruth, cols=starts_with("SRR"), names_to="sample", values_to = "count")

                    do.call(rbind, lapply(c("salmon_out", "kallisto_out", "rsem_out"), function(quant){
                        print(quant)
                        quant_pattern <- ifelse(quant=="salmon_out", "quant.sf", ifelse(quant=="kallisto_out", "abundance.tsv", ".isoforms.results"))
                        txcol <- ifelse(quant=="salmon_out", "Name", ifelse(quant=="kallisto_out", "target_id", "transcript_id"))

                        res_files <- list.files(paste0(sprintf(outdir, seq,rep,R, noise),"/",quant))
                        df_g <- do.call("cbind", lapply(res_files, function(x){
                            quantfiles <- ifelse(quant=="rsem_out", sprintf("%s%s",x,quant_pattern), sprintf("%s", quant_pattern))
                            quantcol <- ifelse(quant=="rsem_out", "expected_count", ifelse(quant=="salmon_out", "NumReads", "est_counts"))

                            tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R, noise), "/", quant, "/", x, "/" ,quantfiles), sep="\t")
                        
                            thistmp <- tmp[,quantcol]
                            thistmp
                        }))
                        df_g <- as.data.frame(df_g)
                        
                        quanttmp <- ifelse(quant=="rsem_out", sprintf("%s/%s%s",res_files[1],res_files[1],quant_pattern), sprintf("%s/%s", res_files[1], quant_pattern))

                        tmp <- read.csv(paste0(sprintf(outdir, seq, rep,R, noise), "/", quant, "/" ,quanttmp), sep="\t")
                        
                        colnames(df_g) <- res_files
                        df_g$tx <- lapply(tmp[,txcol], function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
                        
                        quant_long <- pivot_longer(df_g, cols=starts_with("SRR"), names_to="sample", values_to="count")
                        
                        df_gt <- full_join(groundtruth_long, quant_long, by=c("tx","sample"))
                        df_gt <- df_gt %>% na.omit()
            
                        corr_res <- data.frame('groundtruth'= df_gt$count.x,
                                                'estcount'=df_gt$count.y, 
                                                'quant_tool'=gsub("_out", "", quant),
                                                'seq'=seq,
                                                'rep'=rep,
                                                'sim'=R, 
                                                'tx'=df_gt$tx,
                                                'noise'=noise)
                        
                        rm(df_g)
                        return(corr_res)
                    }))
                }))
            }))
        }))
}))
quantaberr %>% head

png("mae_quant.png", width=900, height=500)
ggplot(quantaberr, aes(x=quant_tool, y=log2(mae)))+
    geom_boxplot(fill="#5cbebb")+
    facet_wrap(rep~seq, labeller=seq_label, scales="free_y")+
    theme_light()+theme(axis.text=element_text(size=20), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        scale_fill_d3(name="Quantification", palette = c("category20"))+
        xlab("Quantification")+
        ylab("RMSE")
dev.off()

png("pairwisetx.png", width=900, height=500)
ggplot(quantaberr, aes(x=groundtruth, y=estcount))+
    geom_point()+
    facet_wrap(quant_tool~seq, labeller=seq_label, scales="free_y")+
    theme_light()+theme(axis.text=element_text(size=20), axis.title = element_text(size=20),
        strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15))+
        scale_fill_d3(name="Quantification", palette = c("category20"))+
        xlab("Quantification")+
        ylab("RMSE")
dev.off()
