library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
path <- args[2]
meta <- args[3]
con1 <- args[4]
con2 <- args[5]
con <- c(con1, con2)
print(con)

jcgene <- fread(paste0(outdir, "/results/junctionseqallGenes.results.txt.gz"))
jcgene <- as.data.frame(jcgene)

resgene <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)))
restx <- read_tsv(paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)))

kresgene <- read_tsv(paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)))
krestx <- read_tsv(paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)))



## gene re
jcgene_g <- data.frame(feature_id=jcgene$geneID, junctionseq=jcgene$geneWisePadj) %>% unique()
resgene <- full_join(resgene, jcgene_g, by="feature_id")
kresgene <- full_join(kresgene, jcgene_g, by="feature_id")
#Function to expand data

## each junction could have two or more transcript + p-values
## 
jcgene_tx <- data.frame(feature_id=jcgene$transcripts, junctionseq=jcgene$padjust_noFilter)

jcgene_tx <- jcgene_tx %>% 
    mutate(feature_id=strsplit(feature_id, ";")) %>% 
    unnest(feature_id) %>% group_by(feature_id) %>% summarise(junctionseq=min(junctionseq))

restx <- full_join(restx, jcgene_tx, by="feature_id")
krestx <- full_join(krestx, jcgene_tx, by="feature_id")

write.table(resgene, paste0(outdir, sprintf("/results/salmon_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
write.table(restx, paste0(outdir, sprintf("/results/salmon_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

write.table(kresgene, paste0(outdir, sprintf("/results/kal_res_gene_%s_%s.txt", con1, con2)), row.names = FALSE, sep="\t")
write.table(krestx, paste0(outdir, sprintf("/results/kal_res_tx_%s_%s.txt", con1, con2)), row.names=FALSE, sep="\t")

