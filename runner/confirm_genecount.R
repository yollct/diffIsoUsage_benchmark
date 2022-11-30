library(tidyverse)
library(countToFPKM)
setwd("/nfs/home/students/chit/is_benchmark")

###
samples_star <- list.files("./alignments", pattern = "sample")


star_output_counts <- list()
count<-1
for (sample in samples_star){
  reads <- read.csv(sprintf("./alignments/%s/genes.fpkm_tracking", sample), sep="\t")
  samplename <- sample
  genecounts <- reads %>% select(gene_id, FPKM)
  genecounts <- genecounts[grepl("ENS", genecounts$gene_id),c(1,2)]
  colnames(genecounts) <- c("gene", "star_gene_count")
  genecounts$sample <- samplename
  
  star_output_counts[[count]] <- genecounts
  count<-count+1
}

all_star <- do.call(rbind, star_output_counts)
all_star$tmp_gene_sam <-paste0(all_star$gene, "_", all_star$sample)


load("./simulated_reads/sim_counts_matrix.rda")
counts_matrix <- as.data.frame(counts_matrix)
counts_matrix$gene <- row.names(counts_matrix)
counts_matrix$gene <- lapply(counts_matrix$gene, function(x){strsplit(x, " ")[[1]][1]}) %>% unlist
counts_matrix$gene <- lapply(counts_matrix$gene, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist

genelength <- lapply(row.names(counts_matrix), function(x){as.numeric(strsplit(strsplit(x, " ")[[1]][3], ":")[[1]][5]) - as.numeric(strsplit(strsplit(x, " ")[[1]][3], ":")[[1]][4])}) %>% unlist

fpkm <- function(x){
  return(x / (genelength * sum(x)) * 1e9)}

fpkmcounts <- as.data.frame(lapply(counts_matrix %>% select(-gene), fpkm))
fpkmcounts$gene <- counts_matrix$gene

groundtruth <- as.data.frame(fpkmcounts)  %>% group_by(gene) %>% dplyr::summarise(across(starts_with("sample"), sum))
groundtruth <- groundtruth %>% gather(key="sample", value="sim_gene_counts", starts_with("sample"))
groundtruth$gene <- lapply(groundtruth$gene, function(x){strsplit(x, " ")[[1]][1]})
groundtruth$tmp_gene_sam <- paste0(groundtruth$gene, "_", groundtruth$sample)

allplot <- inner_join(all_star, groundtruth, by="tmp_gene_sam")

g <- ggplot(allplot, aes(star_gene_count, sim_gene_counts)) +
  geom_point()+
  facet_wrap(.~sample.x)

ggsave(g, file="./results/confirm_genec_star.png")


### for cufflinks
load("./simulated_reads/sim_counts_matrix.rda")
counts_matrix <- as.data.frame(counts_matrix)
counts_matrix$gene <- row.names(counts_matrix)
counts_matrix$gene <- lapply(counts_matrix$gene, function(x){strsplit(x, " ")[[1]][1]}) %>% unlist
counts_matrix$gene <- lapply(counts_matrix$gene, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist

genelength <- lapply(row.names(counts_matrix), function(x){as.numeric(strsplit(strsplit(x, " ")[[1]][3], ":")[[1]][5]) - as.numeric(strsplit(strsplit(x, " ")[[1]][3], ":")[[1]][4])}) %>% unlist
meanfragment <- rep(199, 6)

fpkm <- function(x){
  return(x / (genelength * sum(x)) * 1e6)}

fpkmcounts <- as.data.frame(lapply(counts_matrix %>% select(-gene), fpkm))
fpkmcounts$gene <- counts_matrix$gene
iso_groundtruth <- fpkmcounts %>% gather(key="sample", value="sim_isoform_counts", starts_with("sample"))
iso_groundtruth$gene <- lapply(iso_groundtruth$gene, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
iso_groundtruth$tmp_iso_sam <- paste0(iso_groundtruth$gene, "_", iso_groundtruth$sample)

cufflink_count <- read.csv("./simulated_reads/cufflink_isocounts.csv", sep="\t")
cufflink_count_long <- cufflink_count %>% gather(key="sample", value="isoform_counts", starts_with("sample"))
cufflink_count_long$tmp_iso_sam <- paste0(cufflink_count_long$transcript_id ,"_", cufflink_count_long$sample )

iso_allplot <- inner_join(cufflink_count_long, iso_groundtruth, by=c("tmp_iso_sam"))
g <- ggplot(iso_allplot, aes(isoform_counts, sim_isoform_counts)) +
  geom_point()

ggsave(g, file="./results/confirm_iso_cufflink.png")

### for salmon
load("./simulated_reads/sim_counts_matrix.rda")
counts_matrix <- as.data.frame(counts_matrix)
counts_matrix$gene <- lapply(row.names(counts_matrix), function(x){strsplit(x, " ")[[1]][1]}) %>% unlist
iso_groundtruth <- counts_matrix %>% gather(key="sample", value="sim_isoform_counts", starts_with("sample"))
iso_groundtruth$gene <- lapply(iso_groundtruth$gene, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
iso_groundtruth$tmp_iso_sam <- paste0(iso_groundtruth$gene, "_", iso_groundtruth$sample)


salmon_count <- read.csv("./simulated_reads/salmon_count.csv", sep="\t")
salmon_count_long <- salmon_count %>% gather(key="sample", value="isoform_counts", starts_with("sample"))
salmon_count_long$GENE <- lapply(salmon_count_long$GENE, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist
salmon_count_long$tmp_iso_sam <- paste0(salmon_count_long$GENE ,"_", salmon_count_long$sample )

iso_allplot <- inner_join(salmon_count_long, iso_groundtruth, by=c("tmp_iso_sam"))
g <- ggplot(iso_allplot, aes(isoform_counts, sim_isoform_counts)) +
  geom_point()

ggsave(g, file="./results/confirm_iso_salmon.png")
