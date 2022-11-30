install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",
repos = NULL,
type = "source");
BiocManager::install("https://bioconductor.org/packages/release/bioc/src/contrib/DESeq2_1.14.0.tar.gz")
library(JunctionSeq)
library(DESeq2)
source("/nfs/home/students/chit/is_benchmark/scripts/junctionseq_code.R")

args <- commandArgs(trailingOnly=TRUE)
outdir <- "/MOUNT"
path <- args[2]
meta <- args[3]
index <- args[4]

# outdir <- "/nfs/scratch/chit/simulated_real/single2/"
# path <- "/nfs/home/students/chit/is_benchmark"
# meta <- "/nfs/scratch/chit/simulated_real/single2/meta.txt"
# index <- "/nfs/scratch/chit/ref/ens98_star_rsem"

metadata <- read.csv("/MOUNT/meta.txt", sep="\t")

countFiles <- paste0(outdir, "/qort/", metadata$sample_id,
"/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz")

jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                                sample.names = metadata$sample_id,
                                condition=factor(metadata$group),
                                flat.gff.file = paste0(index, "/JunctionSeq.flat.gff.gz"),
                                nCores = 8,
                                analysis.type = "junctionsAndExons"
                                )

#Generate the size factors and load them into the JunctionSeqCountSet:
jscs <- estimateJunctionSeqSizeFactors(jscs)
jscs <- estimateJunctionSeqDispersions(jscs, nCores = 8)
jscs <- fitJunctionSeqDispersionFunction(jscs)
jscs <- testForDiffUsage(jscs, nCores = 8)
jscs <- estimateEffectSizes( jscs, nCores = 8)

writeCompleteResults(jscs,
outfile.prefix=paste0(outdir, "/results/junctionseq."),
save.jscs = TRUE)

#docker run -v $(pwd):/MOUNT --user $(id -u):$(id -g) --rm --name 'jcseq-run' jcseq --config /nfs/home/students/chit/is_benchmark/single_config.sh
