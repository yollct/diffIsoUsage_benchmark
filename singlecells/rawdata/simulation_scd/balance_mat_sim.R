library(tidyverse)
library(VGAM)
library(Boom)
library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)
library(scDesign3)
library(gridExtra)
library(ggplot2)
library(splatter)

theme_set(theme_bw())
sessionInfo()

# example_data <- readRDS("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/example_data.rds")
rep <- c(50, 100, 200, 500, 700) # how many replicates I want

noise_level <- 0

## specify data paths ##
readfilesdir <- "/nfs/proj/is_benchmark/singlecells/rawdata/smarts2_salmon_est_count.h5ad" # directory to rsem output (as input)
outdir <- paste0("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd")
 # directory to output (output simulated counts)
meta <- paste0(paste0("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd"), "/meta.txt")

meta <- "/nfs/proj/is_benchmark/singlecells/rawdata/" # meta file (tab-separated)
#### this is how the meta file looks like #####
#### first column : sample_id (in your case cell id)
#### second column : group for comparison 
# sample_id	group
# SRR12682098	N
# SRR12682103	N
# SRR12682108	N
# SRR12682113	N
# SRR12682147	T
# SRR12682151	T
# SRR12682155	T
# SRR12682159	T
#############################
unlink(paste0(outdir, "/sim_data_balance"), recursive = TRUE)
dir.create(paste0(outdir, "/sim_data_balance"))

lapply(rep, function(r){
    # if (file.exists(paste0(outdir, sprintf("/sim_data_balance/sim_balance_2ct_rep%s.rds", r)))){
    #     return()
    # }
    print("running")
    print(r)
        ######read single cell data and define clusters
    sce <- readH5AD(readfilesdir)
    geneid <- read.csv("/nfs/proj/is_benchmark/singlecells/rawdata/smarts2_rsem_est_count.h5ad.genetx.txt",sep="\t")
    cellid <- read.csv("/nfs/proj/is_benchmark/singlecells/rawdata/sm2_cellid.out", header=F)
    #celltypes <- 
    pbmc <- CreateSeuratObject(counts=assay(sce,"X"), project="pbmcsm3")
    row.names(pbmc) <- geneid$isoform_id
    colnames(pbmc) <- cellid$V1
    pbmc <- NormalizeData(pbmc)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(pbmc)
    print("here1")
    pbmc <- ScaleData(pbmc, features = all.genes)
    print("here2")
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
    pbmc <- FindNeighbors(pbmc, dims = 1:10)
    pbmc <- FindClusters(pbmc, resolution = 0.5)
    
    #######
    metadata <- data.frame(sample_id=colnames(pbmc),
                            group=Idents(pbmc),
                            rmeans = colMeans(GetAssayData(pbmc[['RNA']], layer="counts")))
    metadata$group <- as.character(metadata$group)
    metadata <- metadata %>% filter(group==0 | group ==1)
    metadata <- metadata %>% arrange(group) %>% group_by(group) %>% sample_n(size = r)
 # sort rows based on group and sample_id column in ascending order
    unigroup <- unique(metadata$group)
    group1n <- lapply(unigroup, function(x){sum(metadata$group==x)}) %>% unlist ### n replicates (number of cells in each clusters) ()
    # Extrahiere die Teile der sample_id nach dem "/"
    # new_values <- sub(".*/", "", metadata$sample_id)

    subpbmc <- subset(pbmc, cells=metadata$sample_id)

    subpbmc <- NormalizeData(subpbmc)
    subpbmc <- FindVariableFeatures(subpbmc, selection.method = "vst", nfeatures = 2000)
    top2000 <- head(VariableFeatures(subpbmc), 2000)

    cellmeans <- rowMeans(GetAssayData(subpbmc, layer="counts"))
    nonzero <- row.names(subpbmc)[cellmeans>5]

    nr <- dim(GetAssayData(subpbmc, layer="counts"))[1]
    nc <- dim(GetAssayData(subpbmc, layer="counts"))[2]

    sce <- as.SingleCellExperiment(subpbmc)
    sce <- sce[nonzero,]

    #write.table(allsrr, file = "allsrr.tsv", sep = "\t", quote = FALSE, col.names = NA)
    allsrr <- GetAssayData(pbmc[['RNA']], layer="counts")
    # Assuming 'allsrr' is your count matrix
    #pseudovalue <- 1

    # Add the pseudovalue to each element
    #allsrr <- allsrr + pseudovalue
    #write.table(allsrr, file = "allsrrpseudo.tsv", sep = "\t", quote = FALSE, col.names = NA)
    ## get binary indicators for mean isoform expression (if mean expression > 0, indicates as 1)
    #exprs <- matrix(rep(t(lapply(rowMeans(allsrr), function(x){ifelse(x>0, 1, 0)} %>% unlist)),2), nrow=nrow(allsrr)) 

    ### read one count data and get the transcript id from real data
    srr <- data.frame(gene_id = geneid$gene_id, transcript_id=geneid$isoform_id, expected_count=rowMeans(allsrr))
    originorder <- row.names(subpbmc)

    srr$rowmeans <- rowMeans(allsrr)
    srr$argsort <- 1:nrow(srr)
    srr <- srr %>% group_by(gene_id) %>% arrange(desc(expected_count), .by_group=T)
    #by_group ensures that sorting is done within each group
    srr <- srr %>% filter(transcript_id %in% nonzero)
    toreorder <- match(nonzero, srr$transcript_id)

    srr <- srr[toreorder,]
    ## get n of transcripts per gene
    metadata <- as.data.frame(metadata)
    row.names(metadata) <- metadata$sample_id
    metadata$group <- factor(metadata$group)

    ## get the gene expression means for simulation
    genecount <- srr %>% group_by(gene_id) %>% mutate(mean_gene_c=sum(rowmeans)) ## mean of expression for genes across all isoforms and samples
    genecount <- genecount %>% mutate(iso_ratio=ifelse(is.na(rowmeans/mean_gene_c), 0, rowmeans/mean_gene_c))
    # This column represents the ratio of expression of each isoform to the mean expression of the gene it belongs to.
    ntrs <- genecount %>% group_by(gene_id) %>% summarise(ntr=n_distinct(transcript_id)) # number of transcripts unique for each gene
    genemeans <- genecount %>% group_by(gene_id) %>% summarise(gm=as.integer(sum(expected_count))) #sum of expression count (means) for each gene

    #### CHANGED TO 20 MIO READS #####
    # genecountcon1 <- rowMeans(data.frame(allsrr[,1:as.integer(rep)]))
    # expectedgenecountcon1 <- genecountcon1 / sum(genecountcon1) * 150000 #change here for Single cell
    # st <- nrow(metadata)-as.integer(rep)+1
    # genecountcon2 <- rowMeans(data.frame(allsrr[,st:nrow(metadata)]))

    rm(example_data)
    rm(example_marginal)
    rm(example_copula)
    ##### scDEsign3

    example_data <- construct_data(
        sce = sce,
        assay_use = "counts",
        celltype = "seurat_clusters",
        pseudotime = NULL,
        spatial = NULL,
        other_covariates = NULL,
        corr_by = "1"
    )

    print("fitting marginal")
    example_marginal <- fit_marginal(
        data = example_data,
        predictor = "gene",
        mu_formula = "seurat_clusters",
        sigma_formula = "1",
        family_use = "nb",
        n_cores = 1,
        usebam = FALSE
    )

    print("fitting copula")
    example_copula <- fit_copula(
        sce = sce,
        assay_use = "counts",
        marginal_list = example_marginal,
        family_use = "nb",
        copula = "gaussian",
        n_cores = 1,
        input_data = example_data$dat
    )
    saveRDS(example_marginal, "example_margin.rds")
    saveRDS(example_data, "example_data.rds")

    print("extracting parameters")
    example_para <- extract_para(
        sce = sce,
        marginal_list = example_marginal,
        n_cores = 1,
        family_use = "nb",
        new_covariate = example_data$newCovariate,
        data = example_data$dat
    )
    print("done")
    # saveRDS(example_para, "example_para.rds")

    ############################################ change mean_mat
    oldmat <- t(example_para$mean_mat)

    oldmat <- data.frame(oldmat)
    oldmat$isoform_id <- row.names(oldmat)
    oldmat <- left_join(oldmat, geneid, by=c("isoform_id"))
    groupings <- oldmat$gene_id
    oldmat <- oldmat %>% select(-isoform_id, -X)

    newmat <- oldmat %>% group_by(gene_id) %>% mutate(across(row.names(example_para$mean_mat),sum)) %>% ungroup() 
    sim_geneid <- newmat$gene_id
    newmat <- newmat %>% dplyr::select(-gene_id)
    
     ##### filter genes###
    zeroinflat <- oldmat[!rowMeans(is.na(t(example_para$zero_mat)))*100 < 50,]$gene_id

    negatives <- unique(c(zeroinflat))

    neg_id <- genemeans$gene_id %in% negatives

    #### CHANGES TO DTU ONLY ####

    ## determine isoform switch or differential isoform usage for each gene
    n <- 1
    geneinfo <- do.call(rbind, lapply(ntrs %>% dplyr::select(ntr) %>% unlist, function(x){
        if (x>1 || n %in% neg_id){
            if (runif(1)>noise_level){ ### runif generate 0-1, if > noise_level (e.g. 10% of noise), simulate if differential events
                if (runif(1)<0.5){
                    return(c(x,T)) #Inside the differential events simulation, it randomly assigns one of three possibilities (isoform switch scenarios) based on additional random numbers (runif(1)). The three possibilities are: (1) switching to the first isoform, (2) switching to the second isoform, or (3) no switching.
                } else {
                        return(c(x,F))
                }
            } else {
                return(c(x,F))
            }
        } else {
            return(c(x,F))
        }
        n<-n+1#The results of each iteration are combined using do.call(rbind, ...) to create a matrix geneinfo where each row represents a gene and columns represent the number of transcripts (ntr) and the presence of differential events for each scenario.
    }))
     ### geneinfo : 
    ## ntr: number of transcripts
    ## indicates of events 
    ## DTU : isoform switch 
    ## DTE : differential transcript expression
    ## DTEDTU : mixture of both
    # "ntr"	"DTE"	"DTU"	"DTEDTU"	"feature_id"
    # 5	1	0	0	"ENSG00000000003"
    # 2	0	1	0	"ENSG00000000005"
    # 6	0	1	0	"ENSG00000000419"
    # 5	1	0	0	"ENSG00000000457"
    # 9	1	0	0	"ENSG00000000460"
    # 7	1	0	0	"ENSG00000000938"
    # 6	1	0	0	"ENSG00000000971"
    # 3	0	1	0	"ENSG00000001036"
    # so delete DTA and dtedtu? because we only need DTU which stands for isoform switch


    #### CHANGE TO DTU ONLY ####
    ## simulate fold change  (determine fold change for each gene (not transcript yet))
    n<-1
    isofc <- c()

    fold_change <- do.call(rbind, lapply(ntrs %>% dplyr::select(ntr) %>% unlist, function(x){
        tmp <- geneinfo[n,2] 
    
        # In summary, the code randomly selects a single value from the sequence (2, 3, 4, 5) and assigns it to the variable fc. The purpose of this code may be to generate a random fold change value for further calculations or simulations.
        fc <- sample(seq(2, 5), size=1)

        if (tmp==0) {
            ### type = 0 : not differentially changing
            isofc <- c(isofc, 1)
            n<<-n+1
        
            c1 = matrix(rep(1,x)) #why rep?
            c2 = matrix(rep(1,x))
            m = cbind(c1,c2)
            return(m)
        }  else if (tmp==1) {
            ### type = 1 : DTU == isoform switch (differentially expression of 2 transcript)
    
            if (x>2){
                m = matrix(c(fc, rep(1, x - 1)))
                c2 = matrix(c(1, fc, rep(1,x-2)))
                m = cbind(m, c2)
                
            } else if (x==2) {
                m = matrix(c(fc ,1))
                c2 = matrix(c(1, fc))
                m = cbind(m, c2)
            } else {
                m = matrix(c(1,1), ncol=2)
            }
            n<<-n+1
            
            # for (i in 2:2) {
            # c2 = sample(m[,1])
            # while (identical(m[,1], c2)) {
            #     c2 <- sample(m[,1])
            # }
            
            isofc <- c(isofc,fc)
            return(m)
        }
    }
    ))
    #### fold change matrix:
    #        c2
    # [1,] 1  4
    # [2,] 1  1
    # [3,] 1  1
    # [4,] 1  1
    # [5,] 1  1
    # [6,] 5  1

    ### simulate isoform ratio using dirichlet distribution

    switch_iso <- c()
    lastgene <- 1
    last_txid <- 0
    matrix_of_sums <- lapply(1:nrow(genemeans), function(x){
        ntr <- ntrs$ntr[x]
        genesum <- genemeans$gm[x]
        tmp <- geneinfo[x,2]

        # dispgene <- dispd[(last_txid+1):(last_txid+ntr),]$disp
        # zero_disp <- which(dispgene==0)

        # switch isoform ratio if the event is isoform switch
        if (ntr>2){
            # if number of transcript > 2
            si <- sample(seq(2, ntr), size=1)
            while (si > ntr){
                si <- sample(seq(2, ntr), size=1)
            }
            alphas <-  c(rep(10, si), rep(1, ntr-si))
            # alphas[zero_disp] <- 0.00000001
            originalgeneratio <- rdirichlet(1, alphas)

            if (si==ntr){
                flip_ori_generatio <- sample(originalgeneratio)
            } else {
            flip_ori_generatio <- c(sample(originalgeneratio[1:si]), originalgeneratio[(si+1):ntr])
            }
        } else {
            si <- ntr
            alphas <- rep(10, ntr)
            # alphas[zero_disp] <- 0.00000001
            originalgeneratio <- rdirichlet(1, alphas)
            flip_ori_generatio <- 1-originalgeneratio
        }
        
        switch_iso <- c(switch_iso, si)
        flip_ori_generatio[flip_ori_generatio==1] <- 0

        lastgene <<- lastgene+ntr
        if (tmp==0) {
            c1<-matrix(originalgeneratio)
            c2<-matrix(originalgeneratio)

            m<-cbind(c1,c2)
            return(m)
        } else if (tmp==1) {
            ### switch isoform ratio if events are DTU or DTEDTU (with isoform switch events)
            
            m = matrix(originalgeneratio)
            c2 = matrix(flip_ori_generatio)
            
            
            m = cbind(m, c2)

            return(m) 
        }

        last_txid <<- ntr
    })

    gene_ratio <- Reduce("rbind", matrix_of_sums)
    #saveRDS(gene_ratio, paste0(readfilesdir, "/gene_ratio_collapsed.rda"))
    ### gene_ratio 
    #           [,1]      [,2]
    # [1,] 0.4543437 0.4543437
    # [2,] 0.2089931 0.2089931
    # [3,] 0.0000000 0.0000000
    # [4,] 0.3366632 0.3366632
    # [5,] 0.0000000 0.0000000
    # [6,] 0.5070024 0.4929976


    ### fold change indicating matrix with replicates
    fold_change_rep <- do.call(cbind, lapply(1:length(group1n), function(x){
        matrix(rep(as.numeric(t(fold_change[,x])), each = group1n[x]), nrow = nrow(fold_change), byrow = TRUE)
    }))
    print(fold_change_rep %>% dim)
    #### output: (4 replicates)
    #      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
    # [1,]    1    1    1    1    4    4    4    4
    # [2,]    1    1    1    1    1    1    1    1
    # [3,]    1    1    1    1    1    1    1    1
    # [4,]    1    1    1    1    1    1    1    1
    # [5,]    1    1    1    1    1    1    1    1
    # [6,]    5    5    5    5    1    1    1    1

    ### isoform ratio indicating matrix with replicates
    gene_ratio_rep <- do.call(cbind, lapply(1:length(group1n), function(x){
        matrix(rep(as.numeric(t(gene_ratio[,x])), each = group1n[x]), nrow = nrow(gene_ratio), byrow = TRUE)
    }))
    print(gene_ratio_rep %>% dim)

    print(oldmat %>% dim)
    #### output 4 replicates
    #           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
    # [1,] 0.4543437 0.4543437 0.4543437 0.4543437 0.4543437 0.4543437 0.4543437
    # [2,] 0.2089931 0.2089931 0.2089931 0.2089931 0.2089931 0.2089931 0.2089931
    # [3,] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    # [4,] 0.3366632 0.3366632 0.3366632 0.3366632 0.3366632 0.3366632 0.3366632
    # [5,] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    # [6,] 0.5070024 0.5070024 0.5070024 0.5070024 0.4929976 0.4929976 0.4929976
    #           [,8]
    # [1,] 0.4543437
    # [2,] 0.2089931
    # [3,] 0.0000000
    # [4,] 0.3366632
    # [5,] 0.0000000
    # [6,] 0.4929976

    ## expected gene count (gene mean )
    # baseexp <- matrix(expectedgenecountcon1)
    # baseexp <- cbind(baseexp, replicate(sum(group1n)-1,baseexp[,1]))
    # baseexp <- as.data.frame(baseexp)
    # baseexp$gene_id <- srr$gene_id
    # baseexp <- baseexp %>% group_by(gene_id) %>% mutate(across(colnames(allsrr), sum)) %>% ungroup %>% dplyr::select(-gene_id)

    ### generate simulated isoform expressions that will be used in negative binomial dis.
    #### the final value for each transcript combines the mean gene expression * isoform ratio * fold change
    final <- gene_ratio_rep * fold_change_rep
    # example_para <- readRDS("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/example_para.rds")
    
    
    print("chcek dim")
    print(newmat %>% dim)
    print(final %>% dim)
    example_para$mean_mat <- t(newmat * final)
    ##########################################


    

    example_para$sigma_mat %>% dim

    print("generating new count")
    example_newcount <- simu_new(
        sce = sce,
        mean_mat = example_para$mean_mat,
        sigma_mat = example_para$sigma_mat,
        zero_mat = example_para$zero_mat,
        quantile_mat = NULL,
        copula_list = example_copula$copula_list,
        n_cores = 1,
        family_use = "nb",
        input_data = example_data$dat,
        new_covariate = example_data$newCovariate,
        important_feature = example_copula$important_feature,
        filtered_gene = example_data$filtered_gene
    )
    print("done")
    saveRDS(example_newcount, "example_newcount.rds")

        #### meta data and groundtruth
    colnames(geneinfo) <- c("ntr", "DTU")
    geneinfo <- as.data.frame(geneinfo)

    isoinfo <- genecount %>% 
                    dplyr::select(transcript_id, gene_id, iso_ratio) 
    isoinfo$isoratioC2 <- gene_ratio[,2]
    geneinfo$feature_id <- unique(isoinfo$gene_id) %>% unlist
    isoinfo$fc <- pmax(fold_change[,1], fold_change[,2])
    isoinfo$rowsum <- rowSums(example_newcount)
    isoinfo$DE.status.bin <- ifelse(isoinfo$rowsum==0, 0, ifelse(fold_change[,1] + fold_change[,2]>2, 1, 0))

    # isoinfo <- isoinfo %>% dplyr::mutate(DE.status.bin=ifelse(transcript_id %in% negatives, 0, DE.status.bin))
    isoinfo$diffiso <- round(abs(gene_ratio[,2]-gene_ratio[,1]),1)

    isoinfo <- isoinfo %>% group_by(gene_id) %>% mutate(new_status=sum(DE.status.bin), n_tx=n()) %>% ungroup
    isoinfo <- isoinfo %>% 
                    dplyr::mutate(gene_group = ifelse(n_tx>9, ">9", ifelse(n_tx>=5, "5-9", ifelse(n_tx>=2, "2-4", "1"))))

    # results_path <- paste0(paste0( 
    # dir.create(results_path, recursive = TRUE)
    # dir.create(simmy_path, recursive = TRUE)
    # dir.create(simnk_path, recursive = TRUE)

    saveRDS(example_newcount, paste0(outdir, sprintf("/sim_data_balance/sim_balance_2ct_rep%s.rds", r)))
    write.table(isoinfo, paste0(outdir, sprintf("/sim_data_balance/sim_balance_2ct_rep%s_isoinfo.csv", r)), sep="\t", row.names = F)
    write.table(geneinfo, paste0(outdir, sprintf("/sim_data_balance/sim_balance_2ct_rep%s_geneinfo.csv", r)), sep="\t", row.names = F)
})
# ### column names : sample name or cell ID 
# row.names(allsrr) <- srr$transcript_id
# colnames(allsrr) <- lapply(resultfiles, function(x){strsplit(x, "/")[[1]][[1]]}) %>% unlist

# ### generate simulated isoform expressions for each condition using dispersion from DESeq2 and real gene count
# simdf <- do.call(cbind, lapply(1:ncol(final), function(xx){
#     tmp1 <- final[,xx]
    
#     ## loop for each transcript and apply negative binomial to generate replicates
#     eachcol <- do.call(c, lapply(1:length(tmp1), function(xxx){ 
#         if (dispd$disp[xxx]==0){
#             return(0)
#         } else {
            
#             size<-1/(dispd$disp[xxx])
            
#             if (is.na(size)){
#                 print(size)
#                 print(tmp1[xxx])
#                 print(dispd$disp[xxx]) ####todo check 0 dispersion
#             }
#             rnbinom(n=1, size=size, mu=tmp1[xxx]) ## 
#         }
#     }))
#     eachcol
# }))
# simdf[is.na(simdf)] <- 0



# #### generate RSEM output for each column using the sample name from the real dataset (or cell ID in single cells)
# lapply(1:sum(group1n), function(x){
#     tmp <- data.frame(transcript_id=srr$transcript_id, gene_id=srr$gene_id, effective_length=srr$effective_length)
#     row.names(tmp) <- tmp$transcript_id
    
#     tmp$expected_count <- simdf[,x]
#     scale_factors<-simdf[,x] / tmp$effective_length
#     scale_factors[is.na(scale_factors) | scale_factors==Inf] <- 0
#     tmp_tpm <- (simdf[,x]/tmp$effective_length)*10^6/sum(scale_factors)
#     tmp_tpm[is.na(tmp_tpm) | tmp_tpm==Inf] <- 0
#     tmp_fpkm <- (simdf[,x]/tmp$effective_length)*10^9/sum(simdf[,x])
#     tmp_fpkm[is.na(tmp_fpkm) | tmp_fpkm == Inf]<-0
#     tmp$TPM <- tmp_tpm

#     tmp$FPKM <- tmp_fpkm 
#     sum <- tmp %>% group_by(gene_id) %>% mutate(sumgene=sum(expected_count))
#     isopct <- round(tmp$expected_count/sum$sumgene*100, 2)
#     isopct[is.na(isopct)]<-0
#     tmp$IsoPct <- isopct
#     tmp <- tmp[originorder,]
#     write.table(tmp, paste0(outdir,"/rsem_sim/", paste0(metadata$sample_id[x],".isoforms.results")), sep="\t", quote=F, row.names=F)
#     print("done")
#     print(metadata$sample_id[x])
# })


sessionInfo()