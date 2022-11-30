library(tidyverse)
library(ggplot2) 
library(UpSetR)
path <- "/nfs/home/students/chit/is_benchmark"
outdir <- "/nfs/scratch/chit/simulated_real/single4"

groundtruthfiles <- list.files(paste0(outdir, "/fastq_sim"), pattern="isoforms.results") ### change to result
groundtruth <- do.call("cbind", lapply(groundtruthfiles, function(x){
    tmp<-read.csv(paste0(outdir,"/fastq_sim/", x), sep="\t")
    tmp$TPM
    }))
colnames(groundtruth) <- lapply(groundtruthfiles, function(x){gsub(".sim.isoforms.results", "", x)})
tmp<-read.csv(paste0(outdir,"/fastq_sim/", groundtruthfiles[1]), sep="\t")
groundtruth <- as.data.frame(groundtruth)
groundtruth$tx <- tmp$transcript_id

salmonfiles <- list.files(paste0(outdir, "/salmon_out"), pattern="SRR")
salmondf <- do.call("cbind", lapply(salmonfiles, function(x){
    tmp<- read.csv(paste0(outdir, "/salmon_out/", x, "/quant.sf"), sep="\t")
    tmp$TPM
}))
colnames(salmondf) <- salmonfiles
tmp<- read.csv(paste0(outdir, "/salmon_out/", salmonfiles[1], "/quant.sf"), sep="\t")
rownames(salmondf) <- lapply(tmp$Name, function(x){strsplit(x, "[.]")[[1]][1]})
salmondf <- as.data.frame(salmondf)
salmondf$tx <- lapply(tmp$Name, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist

kalfiles <- list.files(paste0(outdir, "/kallisto_out"), pattern="SRR")
kaldf <- do.call("cbind", lapply(kalfiles, function(x){
    tmp <- read_tsv(paste0(outdir,"/kallisto_out/",x,"/abundance.tsv"))
    tmp$tpm
}))
colnames(kaldf)<-kalfiles
tmp <- read.csv(paste0(outdir,"/kallisto_out/",kalfiles[1],"/abundance.tsv"), sep="\t")
kaldf <- as.data.frame(kaldf)
kaldf$tx <- lapply(tmp$target_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
kaldf %>% head
sum(!kaldf$tx %in% groundtruth$tx)

rsemfiles <- list.files(paste0(outdir, "/rsem_out"),pattern="SRR")
rsemdf <- do.call("cbind", lapply(rsemfiles, function(x){
    tmp <- read_tsv(paste0(outdir,"/rsem_out/",x,"/",x,".isoforms.results"))
    tmp$TPM
}))
colnames(rsemdf)<-rsemfiles
tmp <- read_tsv(paste0(outdir,"/rsem_out/",rsemfiles[1],"/",rsemfiles[1],".isoforms.results"))
rsemdf <- as.data.frame(rsemdf)
rsemdf$tx <- tmp$transcript_id




##### turn to long 
kaldf_long <- pivot_longer(kaldf, cols=starts_with("SRR"), names_to="sample", values_to = "tpm")
#kaldf_long$id <- paste0(kaldf_long$tx,"_",kaldf_long$sample)

salmondf_long <- pivot_longer(salmondf, cols=starts_with("SRR"), names_to="sample", values_to="tpm")
#salmondf_long$id <- paste0(salmondf_long$tx,"_",salmondf_long$sample)

groundtruth_long <- pivot_longer(groundtruth, cols=starts_with("SRR"), names_to="sample", values_to = "tpm")
#groundtruth_long$id <- paste0(groundtruth_long$tx,"_",groundtruth_long$sample)

rsemdf_long <- pivot_longer(rsemdf, cols=starts_with("SRR"), names_to="sample", values_to="tpm")

kal_gt <- inner_join(groundtruth_long, kaldf_long, c("tx","sample"))
sal_gt <- inner_join(groundtruth_long, salmondf_long, by=c("tx","sample"))
rsem_gt <- inner_join(groundtruth_long, rsemdf_long, by=c("tx", "sample"))
ks_gt <- inner_join(kaldf_long, salmondf_long, by=c("tx", "sample"))
#kal_gt$tpm.x <- as.numeric(kal_gt$tpm.x)
#kal_gt$tpm.y <- as.numeric(kal_gt$tpm.y)

png(paste0(outdir, "/results/kallisto_tx.png"))
ggplot(kal_gt, aes(x=tpm.x, y=tpm.y))+
    geom_point()+
    xlab("Groundtruth TPM")+ylab("Kallisto TPM")+
    theme_classic()+
    geom_text(x=4000, y=30000, label=paste0("R^2=",round(cor(kal_gt$tpm.x, kal_gt$tpm.y, method = "spearman"),2)), size=9)
dev.off()

png(paste0(outdir, "/results/salmon_tx.png"))
ggplot(sal_gt, aes(tpm.x, tpm.y))+
    geom_point()+
    xlab("Groundtruth TPM")+ylab("Salmon TPM")+
    theme_classic()+
    geom_text(x=15000, y=30000, label=paste0("R^2=", round(cor(sal_gt$tpm.x, sal_gt$tpm.y, method = "spearman"),2)), size=9)
dev.off()

png(paste0(outdir, "/results/rsem_tx.png"))
ggplot(rsem_gt, aes(tpm.x, tpm.y))+
    geom_point()+
    xlab("Groundtruth TPM")+ylab("RSEM TPM")+
    theme_classic()+
    geom_text(x=15000, y=30000, label=paste0("R^2=", round(cor(rsem_gt$tpm.x, rsem_gt$tpm.y, method = "spearman"),2)), size=9)
dev.off()

png(paste0(outdir, "/results/sal_kal_tx.png"))
ggplot(ks_gt, aes(tpm.x, tpm.y))+
    geom_point()+
    xlab("Kallisto TPM")+ylab("Salmon TPM")+
    theme_classic()+
    geom_text(x=15000, y=30000, label=paste0("R^2=", round(cor(ks_gt$tpm.x, ks_gt$tpm.y, method = "spearman"),2)), size=9)
dev.off()


