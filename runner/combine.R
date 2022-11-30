library(tidyverse)
path <- "/nfs/home/students/chit/is_benchmark"
source(paste0(path,"/runner/groundtruth.R"))
source(paste0(path,"/runner/event_importance.R"))

dexseq_g <- unname(readRDS(paste0(path, "/results/isa_dexseq_0.rds"))[[1]])
drimseq_g <- unname(readRDS(paste0(path, "/results/isa_drimseq_0.rds"))[[1]])
dtu_g <- unname(readRDS(paste0(path, "/results/dtu_0.rds"))[[1]])

genelist <- list(DEXSeq=unname(dexseq_g[[1]]),
    DRIMSeq=unname(drimseq_g[[1]]),
    DTUrtle=unname(dtu_g[[1]]))
ggvenn(genelist)
ggsave(paste0(path, sprintf("/results/overlap.png")), device="png")

ei0.1 <- filter_eventim(cutoff=0.1, type="gene")

pr <- overall_precision_recall(c(dtu_g), groundtruth, IS_tx)
pr
overall_precision_recall(dtu_g[dtu_g %in% ei0.1$gene], groundtruth, IS_tx)

