library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
dir <- "single_50_4_r2_0.1"

kres <- read.csv(paste0("/nfs/scratch/chit/simulated_real/",dir,"/results/kal_res_gene_N_T.txt"),sep="\t")
sres <- read.csv(paste0("/nfs/scratch/chit/simulated_real/",dir,"/results/salmon_res_gene_N_T.txt"),sep="\t")
rres <- read.csv(paste0("/nfs/scratch/chit/simulated_real/",dir,"/results/rsem_res_gene_N_T.txt"),sep="\t")

kres <- kres %>% select(-dturtle.x) %>% unique() %>% dplyr::rename("dturtle"=dturtle.y)
sres <- sres %>% select(-dturtle.x) %>% unique() %>% dplyr::rename("dturtle"=dturtle.y)
rres <- rres %>% select(-dturtle.x) %>% unique() %>% dplyr::rename("dturtle"=dturtle.y)

# kres <- kres %>% select(-dturtle) %>% unique()
# sres <- sres %>% select(-dturtle) %>% unique()
# rres <- rres %>% select(-dturtle) %>% unique()

write.table(kres, paste0("/nfs/scratch/chit/simulated_real/",dir,"/results/kal_res_gene_N_T.txt"), sep="\t", row.names=F, quote=F)
write.table(sres, paste0("/nfs/scratch/chit/simulated_real/",dir,"/results/salmon_res_gene_N_T.txt"), sep="\t", row.names=F, quote=F)
write.table(rres, paste0("/nfs/scratch/chit/simulated_real/",dir,"/results/rsem_res_gene_N_T.txt"), sep="\t", row.names=F, quote=F)


