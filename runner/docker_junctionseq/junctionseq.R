
library(tidyverse)
library(JunctionSeq)
library(DESeq2)
#source("/nfs/home/students/chit/is_benchmark/scripts/junctionseq_code.R")

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
con1 <- args[2]
con2 <- args[3]
con <- c(con1, con2)
print("comparing:")
print(con)
# outdir <- "/nfs/scratch/chit/simulated_real/single2/"
# path <- "/nfs/home/students/chit/is_benchmark"
# meta <- "/nfs/scratch/chit/simulated_real/single2/meta.txt"
# index <- "/nfs/scratch/chit/ref/ens98_star_rsem"

if (!file.exists(paste0(outdir, sprintf("/results/jseq_%s_%s_/jscs.RData", con1, con2)))){
    metadata <- read.csv("/MOUNT/meta.txt", sep="\t")

    metadata <- metadata %>% dplyr::filter(group %in% con)
    metadata$sample_id <- lapply(metadata$sample_id, function(x){gsub(".sra", "", x)}) %>% unlist()
    countFiles <- paste0(outdir, "/qort/", metadata$sample_id,
    "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz")

    jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                                    sample.names = metadata$sample_id,
                                    condition=factor(metadata$group),
                                    flat.gff.file = paste0("/MOUNT", "/JunctionSeq.flat.gff.gz"),
                                    nCores = 12,
                                    analysis.type = "junctionsAndExons"
                                    )

    #Generate the size factors and load them into the JunctionSeqCountSet:
    jscs <- estimateJunctionSeqSizeFactors(jscs)
    jscs <- estimateJunctionSeqDispersions(jscs)
    jscs <- fitJunctionSeqDispersionFunction(jscs)
    jscs <- testForDiffUsage(jscs)
    jscs <- estimateEffectSizes(jscs)

    if (!file.exists(paste0(outdir, sprintf("/results/jseq_%s_%s_/", con1, con2)))){
        dir.create(paste0(outdir, sprintf("/results/jseq_%s_%s_/", con1, con2)))
    }

    writeCompleteResults(jscs,
    outfile.prefix=paste0(outdir, sprintf("/results/jseq_%s_%s_/", con1, con2)),
    save.jscs = TRUE)

}