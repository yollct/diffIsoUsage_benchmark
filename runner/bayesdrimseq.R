#install.packages('LaplacesDemon')
#install.packages('doMC')
source('/nfs/proj/is_benchmark/runner/BayesDRIMSeq.R')
library('tidyverse')
library('optparse')


option_list <- list(
make_option(c("-i", "--input"), default=NA, type='character', help="the full path and the file name with gene to transcript mapping"),
make_option(c("-d", "--dir"), default=NA, type='character', help="data directory"),
make_option(c("-t", "--tool"), default=NA, type="character", help="quantification tool"))

opt = parse_args(OptionParser(option_list=option_list))

print(opt)
unlink(file.path(opt$dir, "BayesDRIMSeq_tmp"), recursive=T)
#count.txt should have particular format. The first column should comtain gene names from which the transcript counts was calculated
#as far as I understand it should be estimated counts
if (opt$tool=="kal"){
    countDataFrame = read.table(opt$input, header = TRUE,sep="\t")
} else {
    countDataFrame = read.table(opt$input, header = TRUE,sep=",")
}
head(countDataFrame)


genes = countDataFrame$gene_id
#deleting name column as the method itself needs only counts
countDataFrame$gene_id <- NULL
countDataFrame$feature_id <- NULL
reps <- ncol(countDataFrame) / 2
reps
#design matrix - in this case 3 from one condition vs 3 the other condition
grouping = as.factor(c(rep('N',reps), rep('T',reps)))
grouping

#the main function. Output_prefix must be different from run to run
myRes <- laplaceDM(
count_data = countDataFrame,
gene_data = genes, 
grouping=grouping,
min_reads_filter = 10,        
nCores = 8, 
lambdaRate = 0.5,
output_prefix = "tmp")

myRes
sigRes = myRes[myRes$fdrTrust > 0.95,]

resgene <- read_tsv(file.path(opt$dir, sprintf("results/%s_res_gene_N_T.txt", opt$tool)))
res_bdrim <- sigRes %>% dplyr::select(geneNames, fdrTrust) %>% dplyr::rename(feature_id=geneNames, BayesDRIMSeq=fdrTrust)
resgene <- full_join(resgene, res_bdrim, by="feature_id")

write.table(resgene, file.path(opt$dir, sprintf("results/%s_res_gene_N_T.txt", opt$tool)), row.names = FALSE, sep="\t")
# write.table(restx, file.path(opt$dir, sprintf("/%s_res_gene_N_T.txt", opt$tool)), row.names = FALSE, sep="\t")
#Retreive DTU genes
# write.csv(myRes,"bds_rsem.csv")



