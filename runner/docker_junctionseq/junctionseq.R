

.libPaths()

# install.packages("/DESeq2_1.10.1.tar.gz", repos=NULL, type="source", lib="/nfs/home/students/chit/R/x86_64-pc-linux-gnu-library/3.3")
library(dplyr)
library(DESeq2)
library(matrixStats)
library(BiocParallel)
library(foreach)
library(doParallel)
library(statmod)


sourceFolder <- function(folder, recursive = FALSE, ...) 
{ 
    files <- list.files(folder, pattern = "[.][rR]$", 
                        full.names = TRUE, recursive = recursive)
    if (!length(files))
        stop(simpleError(sprintf('No R files in folder "%s"', folder)))
    src <- invisible(lapply(files, source, ...))
    message(sprintf('%s files sourced from folder "%s"', length(src), folder))
}

sourceFolder("/JunctionSeq/R")
library(DESeq2)

#source("/nfs/home/students/chit/is_benchmark/scripts/junctionseq_code.R")

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
con1 <- args[2]
con2 <- args[3]
con <- c(con1, con2)
cores <- as.integer(args[4])

biocpar <- cores
print("13 aug")
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
    metadata$sample_id <- lapply(metadata$sample_id, function(x){as.character(x)}) %>% unlist

    jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                                    sample.names = metadata$sample_id,
                                    condition=factor(metadata$group),
                                    flat.gff.file = paste0("/MOUNT", "/JunctionSeq.flat.gff.gz"),
                                    nCores = biocpar,
                                    analysis.type = "junctionsAndExons"
                                    )

    #Generate the size factors and load them into the JunctionSeqCountSet:
    # jscs <- estimateJunctionSeqSizeFactors(jscs)
    # jscs <- estimateJunctionSeqDispersions(jscs, nCores=biocpar)
    # jscs <- fitJunctionSeqDispersionFunction(jscs)
    # jscs <- testForDiffUsage(jscs, nCores=biocpar)
    # jscs <- estimateEffectSizes(jscs, nCores=biocpar)

    if (!file.exists(paste0(outdir, sprintf("/results/jseq_%s_%s_/", con1, con2)))){
        dir.create(paste0(outdir, sprintf("/results/jseq_%s_%s_/", con1, con2)))
    }

    writeCompleteResults(jscs,
    outfile.prefix=paste0(outdir, sprintf("/results/jseq_%s_%s_/", con1, con2)),
    save.jscs = TRUE)

}
# bams = 1:23
# bams
# lapply(bams, function(x){
#     print("here")
#     print(paste0(outdir, sprintf("/results/jseq_%s_%s_%s/jscs.RData", con1, con2, x)))
#     if (!file.exists(paste0(outdir, sprintf("/results/jseq_%s_%s_%s/jscs.RData", con1, con2, x)))){
#         metadata <- read.csv("/MOUNT/meta.txt", sep="\t")

#         metadata <- metadata %>% dplyr::filter(group %in% con)
#         metadata$sample_id <- lapply(metadata$sample_id, function(xx){gsub(".sra", "", xx)}) %>% unlist()
#         countFiles <- paste0(outdir, sprintf("/qort%s/", x), metadata$sample_id,
#         "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz")
#         print("check")
#         print(countFiles[1])
#         metadata$sample_id <- lapply(metadata$sample_id, function(xx){as.character(xx)}) %>% unlist

#         jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
#                                         sample.names = metadata$sample_id,
#                                         condition=factor(metadata$group),
#                                         flat.gff.file = paste0("/MOUNT", "/JunctionSeq.flat.gff.gz"),
#                                         nCores = biocpar,
#                                         analysis.type = "junctionsAndExons"
#                                         )

#         #Generate the size factors and load them into the JunctionSeqCountSet:
#         # jscs <- estimateJunctionSeqSizeFactors(jscs)
#         # jscs <- estimateJunctionSeqDispersions(jscs, nCores=biocpar)
#         # jscs <- fitJunctionSeqDispersionFunction(jscs)
#         # jscs <- testForDiffUsage(jscs, nCores=biocpar)
#         # jscs <- estimateEffectSizes(jscs, nCores=biocpar)

#         if (!file.exists(paste0(outdir, sprintf("/results/jseq_%s_%s_%s/", con1, con2, x)))){
#             dir.create(paste0(outdir, sprintf("/results/jseq_%s_%s_%s/", con1, con2, x)))
#         }

#         writeCompleteResults(jscs,
#         outfile.prefix=paste0(outdir, sprintf("/results/jseq_%s_%s_%s/", con1, con2, x)),
#         save.jscs = TRUE)

#     }
# })

# if (!file.exists(paste0(outdir, sprintf("/results/jseq_%s_%s_2/jscs.RData", con1, con2)))){
#     metadata <- read.csv("/MOUNT/meta.txt", sep="\t")

#     metadata <- metadata %>% dplyr::filter(group %in% con)
#     metadata$sample_id <- lapply(metadata$sample_id, function(x){gsub(".sra", "", x)}) %>% unlist()
#     countFiles <- paste0(outdir, "/qort2/", metadata$sample_id,
#     "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz")
#     metadata$sample_id <- lapply(metadata$sample_id, function(x){as.character(x)}) %>% unlist

#     jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
#                                     sample.names = metadata$sample_id,
#                                     condition=factor(metadata$group),
#                                     flat.gff.file = paste0("/MOUNT", "/JunctionSeq.flat.gff.gz"),
#                                     nCores = biocpar,
#                                     analysis.type = "junctionsAndExons"
#                                     )

#     #Generate the size factors and load them into the JunctionSeqCountSet:
#     # jscs <- estimateJunctionSeqSizeFactors(jscs)
#     # jscs <- estimateJunctionSeqDispersions(jscs, nCores=biocpar)
#     # jscs <- fitJunctionSeqDispersionFunction(jscs)
#     # jscs <- testForDiffUsage(jscs, nCores=biocpar)
#     # jscs <- estimateEffectSizes(jscs, nCores=biocpar)

#     if (!file.exists(paste0(outdir, sprintf("/results/jseq_%s_%s_2/", con1, con2)))){
#         dir.create(paste0(outdir, sprintf("/results/jseq_%s_%s_2/", con1, con2)))
#     }

#     writeCompleteResults(jscs,
#     outfile.prefix=paste0(outdir, sprintf("/results/jseq_%s_%s_2/", con1, con2)),
#     save.jscs = TRUE)

# }

# if (!file.exists(paste0(outdir, sprintf("/results/jseq_%s_%s_3/jscs.RData", con1, con2)))){
#     metadata <- read.csv("/MOUNT/meta.txt", sep="\t")

#     metadata <- metadata %>% dplyr::filter(group %in% con)
#     metadata$sample_id <- lapply(metadata$sample_id, function(x){gsub(".sra", "", x)}) %>% unlist()
#     countFiles <- paste0(outdir, "/qort3/", metadata$sample_id,
#     "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz")
#     metadata$sample_id <- lapply(metadata$sample_id, function(x){as.character(x)}) %>% unlist

#     jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
#                                     sample.names = metadata$sample_id,
#                                     condition=factor(metadata$group),
#                                     flat.gff.file = paste0("/MOUNT", "/JunctionSeq.flat.gff.gz"),
#                                     nCores = biocpar,
#                                     analysis.type = "junctionsAndExons"
#                                     )

#     #Generate the size factors and load them into the JunctionSeqCountSet:
#     # jscs <- estimateJunctionSeqSizeFactors(jscs)
#     # jscs <- estimateJunctionSeqDispersions(jscs, nCores=biocpar)
#     # jscs <- fitJunctionSeqDispersionFunction(jscs)
#     # jscs <- testForDiffUsage(jscs, nCores=biocpar)
#     # jscs <- estimateEffectSizes(jscs, nCores=biocpar)

#     if (!file.exists(paste0(outdir, sprintf("/results/jseq_%s_%s_3/", con1, con2)))){
#         dir.create(paste0(outdir, sprintf("/results/jseq_%s_%s_3/", con1, con2)))
#     }

#     writeCompleteResults(jscs,
#     outfile.prefix=paste0(outdir, sprintf("/results/jseq_%s_%s_3/", con1, con2)),
#     save.jscs = TRUE)

# }