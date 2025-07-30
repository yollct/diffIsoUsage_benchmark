library(BANDITS)
library(tximport)
library(optparse)
library(tidyverse)

option_list <- list(
make_option(c("-g", "--gene_to_transcript"), default=NA, type='character', help="the full path and the file name with gene to transcript mapping"),
make_option(c("-d","--dir"),default=NA, type='character', help="the directory with the salmon results folder"),
make_option(c("-s","--salmon"),default=NA, type='character', help="the name of the salmon files directory"),
make_option(c("-c","--ncores"),default=2, type='integer', help="number of cores"),
make_option(c("-e","--seed"),default=61217, type='integer', help="seed"),
make_option(c("-r","--gene_results"),default=61217, type='character', help="where to save transcript results in a csv format"),
make_option(c("-t","--tr_results"),default=61217, type='character', help="where to save transcript results in a csv format"),
make_option(c("-i","--image"),default=NA, type='character', help="where to save image with the name"))

opt = parse_args(OptionParser(option_list=option_list))
print(opt)

#gene to transcript mapping. Should include all transcripts used in the salmon mapping. Format: gene_id transcript_id
gene_tr_id <- read.table(opt$gene_to_transcript,sep="\t",header=T)

#path to the folder with the folder with salmon results
data_dir = opt$dir

#the name of the salmon results folder
salmon_dir = opt$salmon

#create paths to quant.sf files
sample_names = list.files(file.path(data_dir, salmon_dir))
quant_files = file.path(data_dir, salmon_dir, sample_names, "quant.sf")
quant_files
#check that the files exist. Should be all TRUE
file.exists(quant_files)

#Read salmon quant.sf files
txi = tximport(files = quant_files, type = "salmon", txOut = TRUE)
counts = txi$counts

#Create a design matrix
samples_design = data.frame(sample_id = sample_names, group = c("SHAM", "SHAM", "SHAM", "SHAM", "TAC", "TAC", "TAC", "TAC"))
samples_design

#Compute effective length
eff_len = eff_len_compute(x_eff_len = txi$length)
head(eff_len)

#Optional but recommended step to filter out transcripts
transcripts_to_keep = filter_transcripts(gene_to_transcript = gene_tr_id,
                                         transcript_counts = counts, 
                                         min_transcript_proportion = 0.01,
                                         min_transcript_counts = 10, 
                                         min_gene_counts = 20)

#Create path to equivalence classes file and check that they exist. Must be TRUE for all 
#The eq_classes.txt must be gunzipped first
sample_names
equiv_classes_files = file.path(data_dir, salmon_dir, sample_names, "aux_info", "eq_classes.txt")
file.exists(equiv_classes_files)

#Reading all input info (this step might take a while - more than 2 hours for 8 samples with 10 cores)
print('Reading the data')
input_data = create_data(salmon_or_kallisto = "salmon",
                         gene_to_transcript = gene_tr_id,
                         salmon_path_to_eq_classes = equiv_classes_files,
                         eff_len = eff_len, 
                         n_cores = opt$ncores,
                         transcripts_to_keep = transcripts_to_keep)

save.image(opt$image)

#Filter low-expressed genes
print('Filtering genes')
input_data = filter_genes(input_data, min_counts_per_gene = 20)

set.seed(61217)
#The step is optional but recommended
print('Calculation prior precision')
precision = prior_precision(gene_to_transcript = gene_tr_id,
                            transcript_counts = counts, n_cores = opt$ncores,
                            transcripts_to_keep = transcripts_to_keep)


set.seed(61217)
#DTU test with the parameters from the manual
print('DTU testing')
results = test_DTU(BANDITS_data = input_data,
                   precision = precision$prior,
                   samples_design = samples_design,
                   group_col_name = "group",
                   R = 10^4, burn_in = 2*10^3, n_cores = opt$ncores,
                   gene_to_transcript = gene_tr_id)

save.image(opt$image)
results

#Saving results
resgene <- read_tsv(file.path(data_dir, "/results/salmon_res_gene_N_T.txt"))
gene_results = top_genes(results, n = Inf, sort_by_g = "p.value")
ban_gene <- gene_results %>% dplyr::select(Gene_id, adj.p.values) %>% dplyr::rename(BANDITS=adj.p.values, feature_id=Gene_id)
ban_gene$feature_id <- lapply(ban_gene$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist

print('Saving gene results')
resgene <- full_join(resgene, ban_gene, by="feature_id")
# write.csv(gene_results, opt$gene_results)

restx <- read_tsv(file.path(data_dir, "results/salmon_res_tx_N_T.txt"))
transcript_results = top_transcripts(results, n = Inf)
print('Saving transcripts results')
ban_tx <- transcript_results %>% dplyr::select(Transcript_id, adj.p.values) %>% dplyr::rename(feature_id=Transcript_id, BANDITS=adj.p.values)
ban_tx$feature_id <- lapply(ban_tx$feature_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist

restx <- full_join(restx, ban_tx, by="feature_id")
# write.csv(transcript_results, opt$transcript_results)

write.table(resgene,file.path(data_dir, "results/salmon_res_gene_N_T.txt"), row.names = FALSE, sep="\t")
write.table(restx,file.path(data_dir, "results/salmon_res_tx_N_T.txt"), row.names = FALSE, sep="\t")
save.image(opt$image)
