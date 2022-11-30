suppressMessages(library(polyester))
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(rtracklayer)) 
cl <- makeCluster(4)
registerDoParallel(cl)
clusterCall(cl, function() library(magrittr))
path <- "/nfs/home/students/chit/is_benchmark"

args = commandArgs(trailingOnly=TRUE)
n_gene <- as.numeric(args[1])
depth <- as.numeric(args[2])
switches <- args[3]
rep <- as.numeric(args[4])
outdir <- args[5]
# FASTA annotatio
print(sprintf("number of rep : %s", rep))


fasta <- readDNAStringSet('/nfs/data/covid_hscell_4tp/ensembl_99/Homo_sapiens.GRCh38.cdna.all.fa.gz')

#fasta
fasta_gene <- data.frame(transcript_id = lapply(names(fasta), function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist,
                         gene_id = lapply(names(fasta), function(x){strsplit(strsplit(x, " ")[[1]][4], ":")[[1]][2]}) %>% unlist)

fasta_gene_selected <- fasta_gene$gene_id %in% sample(fasta_gene$gene_id, n_gene)
selected_fasta <- fasta[fasta_gene_selected]
writeXStringSet(selected_fasta, paste0(outdir,"/simulated_reads/selected_fasta.fa"))

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM

readspertx = round(depth * width(selected_fasta) / 100)


ntrs <- fasta_gene[fasta_gene_selected,] %>% group_by(gene_id) %>% summarise(ntr = n_distinct(transcript_id)) %>% select(ntr) %>% unlist
fold_change <- do.call(rbind, lapply(ntrs, function(x){
  if (x==1){
    return(matrix(rep(1,1), ncol=2))
  } else {
    if (runif(1) < 0.1)
      return(matrix(rep(1,x*2), ncol=2))
    
    fc <- sample(seq(2, 10), size=1)
    
    m = matrix(c(fc, rep(1, x - 1)))
    for (i in 2:2) {
      c2 = sample(m[,1])
      while (identical(m[,1], c2)) {
        c2 <- sample(m[,1])
      }
      m = cbind(m, c2)
    }
    return(m)
  }
}))

set.seed(4)
print("simulating reads...")
simulate_experiment(paste0(outdir,"/simulated_reads/selected_fasta.fa"),
                    reads_per_transcript=readspertx, 
                    num_reps=c(rep,rep), fold_changes=fold_change, outdir=paste0(outdir,'/simulated_reads')) 


library(rtracklayer)
gtf=rtracklayer::import('/nfs/data/covid_hscell_4tp/ensembl_106/Homo_sapiens.GRCh38.106.gtf')
gtf$transcript_v <- paste0(gtf$transcript_id, ".", gtf$transcript_version)
salmoncnt[row.names(salmoncnt) %in% unique(gtf$transcript_v),] %>% dim
 
tt <- lapply(names(fasta), function(x){strsplit(x, "[ ]")[[1]][1]}) %>% unlist

files <- list.files(paste0(outdir,'/simulated_reads'))
for (f in files){
  if (grepl(".fasta", f)){
    file.rename(paste0(outdir, "/simulated_reads/", f.))
  }
}