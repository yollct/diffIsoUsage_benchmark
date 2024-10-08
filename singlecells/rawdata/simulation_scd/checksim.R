library(tidyverse)
library(splatter)
library(zellkonverter)

readfilesdir <- "/nfs/home/students/chit/is_benchmark/singlecells/rawdata/smarts2_salmon_est_count.h5ad" 
rawsce <- readH5AD(readfilesdir)
sce1 <- SingleCellExperiment(assays=list(counts=assay(rawsce,"X")))
logcounts(sce1) <- log1p(counts(sce1))

sim2 <- readRDS("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/simulation_scd/sim_data_balance/sim_balance_2ct_rep50.rds")
sce2 <- SingleCellExperiment(assays = list(counts = sim2))
logcounts(sce2) <- log1p(counts(sce2))





fn <- readRDS("listfalseneg.rds")
fp <- readRDS("listfalsepos.rds")
mets <- c("DTUrtle", "satuRn")
reps <- c('50', '100', '200')


fn_sce <- lapply(reps, function(r){
    lapply(mets, function(m){
        sim2 <- readRDS(sprintf("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/simulation_scd/sim_data_balance/sim_balance_2ct_rep%s.rds", r))
        sce2 <- SingleCellExperiment(assays = list(counts = sim2))  
        eachfn <- fn[[r]][[m]]
        print(eachfn %in% row.names(sce2))
        sce2[eachfn,]
        
    })
})
fn_sce <- unlist(fn_sce, recursive = FALSE)
names(fn_sce) <- Map(paste0, mets, base::rep(reps, each=2), "FN", sep="")

fp_sce <- lapply(reps, function(r){
    lapply(mets, function(m){
        sim2 <- readRDS(sprintf("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/simulation_scd/sim_data_balance/sim_balance_2ct_rep%s.rds", r))
        sce2 <- SingleCellExperiment(assays = list(counts = sim2))  
        eachfp <- fp[[r]][[m]]

        sce2[eachfp,]
        
    })
})
fp_sce <- unlist(fp_sce, recursive = FALSE)
names(fp_sce) <- Map(paste0, mets, base::rep(reps, each=2), "FP", sep="")

tp <- readRDS("testpos_genes.rds")
testp_sce <- lapply(reps, function(r){
    lapply(mets, function(m){
        sim2 <- readRDS(sprintf("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/sim_data_balance1_0/sim_balance_2ct_rep%s.rds", r))
        sce2 <- SingleCellExperiment(assays = list(counts = sim2))  
        eachfn <- tp[[paste0(m,r)]]
        print(paste0(m,r))
        sce2[eachfn,]
        
    })
})
testp_sce <- unlist(testp_sce, recursive = FALSE)
names(testp_sce) <- Map(paste0, mets, base::rep(reps, each=2), "TP", sep="")

comp <- compareSCEs(c(fn_sce, testp_sce))
png("simmeans.png")
comp$Plots$Means
dev.off()

png("simvar.png")
comp$Plots$Variance
dev.off()

png("simmeanvar.png")
comp$Plots$MeanVar
dev.off()

png("ZerosGene.png")
comp$Plots$ZerosGene
dev.off()

png("ZerosCell.png")
comp$Plots$ZerosCell
dev.off()


library(DuoClustering2018)
Zhengmix4eq_sce <- get("sce_filteredExpr10_Zhengmix4eq")(metadata = FALSE)
ngene <- 200
logcounts(Zhengmix4eq_sce) <- log1p(counts(Zhengmix4eq_sce))
zheng_sce <- modelGeneVar(Zhengmix4eq_sce)
chosen <- getTopHVGs(zheng_sce, n = ngene)
example_sce <- Zhengmix4eq_sce
selected_cells <- which(colData(example_sce)$phenoid %in% c("b.cells","regulatory.t"))
example_sce <- example_sce[,selected_cells]
colData(example_sce)$cell_type <- as.factor(colData(example_sce)$phenoid)
example_sce
head(colData(example_sce))
set.seed(123)
example_data <- construct_data(
    sce = example_sce,
    assay_use = "counts",
    celltype = "cell_type",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = NULL,
    corr_by = "1"
  )
example_marginal <- fit_marginal(
    data = example_data,
    predictor = "gene",
    mu_formula = "cell_type",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 2,
    usebam = FALSE
  )
set.seed(123)
example_copula <- fit_copula(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = example_marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = 2,
    input_data = example_data$dat
  )
example_para <- extract_para(
    sce = example_sce,
    marginal_list = example_marginal,
    n_cores = 2,
    family_use = "nb",
    new_covariate = example_data$newCovariate,
    data = example_data$dat
  )
library(tidyverse)
example_para$zero_mat %>% head
rowMeans(t(example_para$mean_mat)) %>% head
rowMeans(is.na(t(example_para$zero_mat))) *100



library(DTUrtle)
library(satuRn)
library(tidyverse)
library(biomaRt)
library(satuRn)
library(SummarizedExperiment)
library(tximport)
library(edgeR)
source('/nfs/proj/is_benchmark/singlecells/testDTU.R')

biocpar <- BiocParallel::MulticoreParam(4)


outdir <- "/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/"

reps <- c(2,3,4,5,6,7) #, 100, 200, 500)#, 500,800)
allfolders <- list.files(outdir, pattern="sim_data_unbalance")

human <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=107)


resgene <- do.call(c, lapply(allfolders, function(thisfolder){
    lapply(reps, function(r){
    # Assuming count_matrix is your single-cell count matrix
    # Replace this with your actual matri
    print("run")
    print(r)
    new_count <- readRDS(sprintf("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/simulation_scd/%s/sim_unbalance_2ct_rep%s.rds", thisfolder, r))

    metadata <- read.csv(sprintf("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/simulation_scd/%s/sim_unbalance_2ct_rep%s_metadata.csv", thisfolder, r), sep="\t")
    metadata$group <- paste0("celltype", metadata$group)
    metadata$folder <- thisfolder
    metadata$reps <- r 
    metadata
})}))

resgene <- Reduce(rbind, resgene)
cellno <- resgene %>% group_by(folder, reps) %>% summarise(count=n())
png("cell_num.png")
ggplot(cellno, aes(x=reps, y=count))+
  geom_histogram(stat="identity")
dev.off()

outdir <- "/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/"

reps <- c(50, 100, 200, 500, 700)#, 100, 200, 500)#, 500,800)
allfolders <- list.files(outdir, pattern="sim_data_balance")
# allfolders <- allfolders[grepl("0$",allfolders)] ### only noise = 0

truepos <- lapply(allfolders, function(thisfolder){
    lapply(reps, FUN=function(r){
    # Assuming count_matrix is your single-cell count matrix
    # Replace this with your actual matri
    print("run")
    print(r)
    new_geneinfo <- read.csv(sprintf("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/%s/sim_balance_2ct_rep%s_geneinfo.csv", thisfolder, r), sep="\t")
    
    truepos <- new_geneinfo[new_geneinfo$DTU==1,]$feature_id 

    # truepos <- new_geneinfo[new_geneinfo$DE.status.bin==1,]$transcript_id
    unique(truepos)
})})
truepos <- unlist(truepos, recursive=F)
names(truepos) <- Map(paste0, base::rep(allfolders, each=length(reps)), "_",reps)

simdata <- data.frame(ntruepos = lapply(truepos, function(x){length(x)}) %>% unlist)
simdata


unbalance_final <- read.csv("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/singlecell_unbalance_finallong.csv", sep=" ")
balance_final <- read.csv("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/singlecell_balance_finallong.csv", sep=" ")
unba_meta <- read.csv("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/metacell_balance_finallong.csv", sep=" ")
ba_meta <- read.csv("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/metacell_balance_finallong.csv", sep=" ")


unbalance_final$type <- "single-cell"
balance_final$type <- "single-cell"
unba_meta$type <- "meta-cell"
ba_meta$type <- "meta-cell"

unba_final <- rbind(unbalance_final, unba_meta)
unbalance_final %>% head
unba_meta %>% head

png("singlecell_0_unbalance_prere.png", res=300, width=3000, height=1000)
ggplot(final_long, aes(x=reps, y=mean, color=tool))+
    geom_point()+
    geom_line()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                 position=position_dodge(0.05))+
    facet_grid(.~measure)+
    scale_color_manual(values=c('#999999','#E69F00'))+
    theme_bw()+xlab("Number of cell types")+ylab("Precision/recall")
dev.off()