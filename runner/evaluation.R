library(tidyverse)
library(ggvenn)
path <- "/nfs/home/students/chit/is_benchmark"
source(paste0(path,"/runner/event_importance.R"))

args <- commandArgs(trailingOnly=TRUE)
foldch <- args[1]
depth <- args[2]

all_pr_files <- list.files(paste0(path,"/results/"), pattern = "pr_df")

allggfiles <- all_pr_files[grepl("gg", all_pr_files) & !grepl("tx", all_pr_files)]
allggpr<- do.call("rbind", lapply(allggfiles, function(x){read.csv(paste0(path,sprintf("/results/%s", x)), sep=" ")}))

allprfiles <- all_pr_files[!grepl("gg", all_pr_files) & !grepl("tx", all_pr_files)]
allpr<- do.call("rbind", lapply(allprfiles, function(x){read.csv(paste0(path,sprintf("/results/%s", x)), sep=" ")}))

alltxprfiles <- all_pr_files[!grepl("gg", all_pr_files) & grepl("tx", all_pr_files)]
alltxpr<- do.call("rbind", lapply(alltxprfiles, function(x){read.csv(paste0(path,sprintf("/results/%s", x)), sep=" ")}))

alltxggprfiles <- all_pr_files[grepl("gg", all_pr_files) & grepl("tx", all_pr_files)]
alltxggpr<- do.call("rbind", lapply(alltxggprfiles, function(x){read.csv(paste0(path,sprintf("/results/%s", x)), sep=" ")}))



allpr$foldchange <- foldch
allpr$seqdepth <- depth

allggpr$foldchange <- foldch
allggpr$seqdepth <- depth

alltxpr$foldchange <- foldch
alltxpr$seqdepth <- depth

alltxggpr$foldchange <- foldch
alltxggpr$seqdepth <- depth

write.table(allpr, paste0(path, sprintf("/results/allpr_%s_%s.csv", foldch, depth)), row.names=F)
write.table(allggpr, paste0(path, sprintf("/results/allggpr_%s_%s.csv", foldch, depth)), row.names=F)
write.table(alltxpr, paste0(path, sprintf("/results/alltxpr_%s_%s.csv", foldch, depth)), row.names=F)
write.table(alltxggpr, paste0(path, sprintf("/results/alltxggpr_%s_%s.csv", foldch, depth)), row.names=F)


dexseq_g <- readRDS(paste0(path, "/results/isa_dexseq_0.rds"))
drimseq_g <- readRDS(paste0(path, "/results/isa_drimseq_0.rds"))
dtu_g <- readRDS(paste0(path, "/results/dtu_0.rds"))

genelist <- list(DEXSeq=unname(dexseq_g[[1]]),
    DRIMSeq=unname(drimseq_g[[1]]),
    DTUrtle=unname(dtu_g[[1]]))
ggvenn(genelist)
ggsave(paste0(path, sprintf("/results/overlap.png")), device="png")
