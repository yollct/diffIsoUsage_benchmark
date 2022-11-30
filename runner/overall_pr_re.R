library(tidyverse)

path<-"/nfs/home/students/chit/is_benchmark"

allfiles <- list.files(paste0(path,"/results"),pattern="allpr")
csvfiles <- allfiles[endsWith(allfiles, "csv")]

allpr <- data.frame()
for (f in csvfiles){
    tmp <- read.csv(paste0(path,"/results/",f), sep=" ")
    allpr <- rbind(allpr, tmp)
}
allpr <- allpr %>% dplyr::filter(tool %in% c("dturtle_0","ISA (DEXSeq)_0", "ISA (DRIMSeq)_0","Cuffdiff", "iso-KTSP"))

allpr$cutoff <- as.character(allpr$cutoff)

g <- ggplot(allpr %>% dplyr::filter(eventim=="no"), aes(recall, precision, color=tool, shape=cutoff))+
    geom_line()+
    geom_point()+
    facet_grid(foldchange~seqdepth)

ggsave(file=paste0(path, sprintf("/results/overall_pr.png")))

g <- ggplot(allpr %>% dplyr::filter(eventim=="yes"), aes(recall, precision, color=tool, shape=cutoff))+
    geom_line()+
    geom_point()+
    facet_grid(foldchange~seqdepth)

ggsave(file=paste0(path, sprintf("/results/overall_pr_eifilter.png")))

allggfiles <- list.files(paste0(path,"/results"),pattern="allggpr")
csvggfiles <- allggfiles[endsWith(allggfiles, "csv")]

allggpr <- data.frame()
for (f in csvggfiles){
    tmp <- read.csv(paste0(path,"/results/",f), sep=" ")
    allggpr <- rbind(allggpr, tmp)
}
######### gene group overall gene level ##############
allggpr <- allggpr %>% dplyr::filter(tool %in% c("dturtle_0","ISA (DEXSeq)_0", "ISA (DRIMSeq)_0", "Cuffdiff", "iso-KTSP"))
allggpr$gene_group <- factor(allggpr$gene_group, levels=c("2-4","5-9",">9"))
g<-ggplot(allggpr %>% dplyr::filter(eventim=="no"), aes(gene_group, recall, color=tool))+
    geom_boxplot(position="dodge")+
    facet_grid(foldchange~seqdepth)+
    ggtitle("Recall (gene level)")
ggsave(file=paste0(path, sprintf("/results/overall_gg_recall.png")))

g<-ggplot(allggpr %>% dplyr::filter(eventim=="no"), aes(gene_group, precision, color=tool))+
    geom_boxplot(position="dodge")+
    facet_grid(foldchange~seqdepth)+
    ggtitle("Precision (gene level)")
ggsave(file=paste0(path, sprintf("/results/overall_gg_precision.png")))

## event importnace 

g<-ggplot(allggpr %>% dplyr::filter(eventim=="yes"), aes(gene_group, recall, color=tool))+
    geom_boxplot(position="dodge")+
    facet_grid(foldchange~seqdepth)+
    ggtitle("Recall (gene level)")
ggsave(file=paste0(path, sprintf("/results/overall_gg_recall_eifilter.png")))

g<-ggplot(allggpr %>% dplyr::filter(eventim=="yes"), aes(gene_group, precision, color=tool))+
    geom_boxplot(position="dodge")+
    facet_grid(foldchange~seqdepth)+
    ggtitle("Precision (gene level)")
ggsave(file=paste0(path, sprintf("/results/overall_gg_precision_eifilter.png")))

############# transcript level overall ##################3
alltxfiles <- list.files(paste0(path,"/results"),pattern="alltxpr")
csvtxprfiles <- alltxfiles[endsWith(alltxfiles, "csv")]

alltxpr <- data.frame()
for (f in csvtxprfiles){
    tmp <- read.csv(paste0(path,"/results/",f), sep=" ")
    alltxpr <- rbind(alltxpr, tmp)
}

############## transcript level gene group ###############
alltxpr$postfilter <- lapply(alltxpr$tool, function(x){strsplit(x, "_")[[1]][2]}) %>% unlist 
#alltxpr$postfilter <- as.numeric(alltxpr$postfilter)
alltxpr$tool <- lapply(alltxpr$tool, function(x){strsplit(x, "_")[[1]][1]}) %>% unlist
alltxpr$cutoff <- as.character(alltxpr$cutoff)


g<-ggplot(alltxpr %>% dplyr::filter(postfilter==0 & eventim=="no"), aes(recall, precision, color=tool, shape=cutoff)) +
    geom_point()+
    facet_grid(seqdepth~foldchange)
ggsave(file=paste0(path, sprintf("/results/overall_tx_pr.png")))

g<-ggplot(alltxpr %>% dplyr::filter(postfilter==0 & eventim=="yes"), aes(recall, precision, color=tool, shape=cutoff)) +
    geom_point()+
    facet_grid(seqdepth~foldchange)
ggsave(file=paste0(path, sprintf("/results/overall_tx_pr_eifilter.png")))

# g<-ggplot(alltxpr %>% dplyr::filter(seqdepth==60 & foldchange=="=4" & cutoff=="0.05"), aes(recall, precision, color=tool, shape=cutoff)) +
#     geom_line()+geom_point()+
#     facet_grid(foldchange~gene_group)
# ggsave(file=paste0(path, sprintf("/results/overall_tx_recall.png")))



###post hoc 
posthoc <- alltxpr %>% dplyr::filter(cutoff==0.05)
posthoc$postfilter <- as.numeric(posthoc$postfilter)
g<-ggplot(posthoc %>% dplyr::filter(eventim=="no"), aes(postfilter, precision, shape=tool, color=tool))+
    geom_point(alpha=0.7)+geom_line()+
    facet_grid(seqdepth~foldchange)
ggsave(file=paste0(path, sprintf("/results/posthoc_precision.png")))

g<-ggplot(posthoc %>% dplyr::filter(eventim=="no"), aes(postfilter, recall, shape=tool, color=tool))+
    geom_point(alpha=0.7)+geom_line()+
    facet_grid(seqdepth~foldchange)
ggsave(file=paste0(path, sprintf("/results/posthoc_recall.png")))

g<-ggplot(posthoc %>% dplyr::filter(eventim=="yes"), aes(postfilter, precision, shape=tool, color=tool))+
    geom_point(alpha=0.7)+geom_line()+
    facet_grid(seqdepth~foldchange)
ggsave(file=paste0(path, sprintf("/results/posthoc_precision_eifilter.png")))

g<-ggplot(posthoc %>% dplyr::filter(eventim=="yes"), aes(postfilter, recall, shape=tool, color=tool))+
    geom_point(alpha=0.7)+geom_line()+
    facet_grid(seqdepth~foldchange)
ggsave(file=paste0(path, sprintf("/results/posthoc_recall_eifilter.png")))
