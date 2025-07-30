library(tidyverse)
library(VGAM)
library(DESeq2)
library(Boom)

args <- commandArgs(trailingOnly=TRUE)
newdir <- args[1]
rep <- args[2]
noise_level <- as.numeric(args[3])


## specify data paths ##
readfilesdir <- "/nfs/scratch/chit/covid_ts/rsem_out/" # directory to rsem output (as input)
outdir <- sprintf("/nfs/scratch/chit/new_simulations/%s/", newdir) # directory to output (output simulated counts)
meta <- "/nfs/scratch/chit/new_simulations/single_temp/meta.txt" # meta file (tab-separated)
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

metadata <- read_tsv(meta)
metadata <- metadata %>% arrange(group, sample_id)
unigroup <- unique(metadata$group)
group1n <- c(as.integer(rep), as.integer(rep)) ### n replicates (number of cells in each group) ()

### list of the relative path of rsem counts (each sample has a sub-directory)
resultfiles <- paste0(metadata$sample_id,"/",metadata$sample_id, ".isoforms.results")

### reading all rsem counts 
allsrr <- do.call(cbind, lapply(resultfiles, function(x){
    tmp <- read_tsv(paste0(readfilesdir, x))
    tmp$expected_count
}))

## get binary indicators for mean isoform expression (if mean expression > 0, indicates as 1)
exprs <- matrix(rep(t(lapply(rowMeans(allsrr), function(x){ifelse(x>0, 1, 0)} %>% unlist)),2), nrow=nrow(allsrr)) 

### read one count data and get the transcript id from real data
srr <- read_tsv(paste0(readfilesdir, resultfiles[1]))
originorder <- srr$transcript_id

### getting the major isoforms in each genes (sorting the counts within each genes)
simmat <- matrix(rep(t(srr$expected_count), each=sum(group1n)), nrow=nrow(srr), byrow=T)

srr$rowmeans <- rowMeans(allsrr)
srr$argsort <- 1:nrow(srr)
srr <- srr %>% group_by(gene_id) %>% arrange(desc(expected_count), .by_group=T)
## get n of transcripts per gene

metadata <- as.data.frame(metadata)
row.names(metadata) <- metadata$sample_id
metadata$group <- factor(metadata$group)
mode(allsrr) <- "integer"

## use DESeq2 to estimate dispersion
dds <- DESeqDataSetFromMatrix(countData = allsrr[srr$argsort,],
                        colData = metadata,
                        design= ~group)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

dispd <- as.data.frame(mcols(dds)) %>% dplyr::select(dispGeneEst) %>% mutate(disp=ifelse(is.na(dispGeneEst), 0, dispGeneEst))

## 
genecount <- srr %>% group_by(gene_id) %>% mutate(mean_gene_c=sum(rowmeans)) ## mean of expression for genes across all isoforms and samples
genecount <- genecount %>% mutate(iso_ratio=ifelse(is.na(rowmeans/mean_gene_c), 0, rowmeans/mean_gene_c))
ntrs <- genecount %>% group_by(gene_id) %>% summarise(ntr=n_distinct(transcript_id))
genemeans <- genecount %>% group_by(gene_id) %>% summarise(gm=as.integer(sum(expected_count)))

genecountcon1 <- rowMeans(data.frame(allsrr[,1:as.integer(rep)]))
expectedgenecountcon1 <- genecountcon1 / sum(genecountcon1) * 50000000
st <- nrow(metadata)-as.integer(rep)+1
genecountcon2 <- rowMeans(data.frame(allsrr[,st:nrow(metadata)]))
expectedgenecountcon2 <- genecountcon2 / sum(genecountcon2) * 50000000

# genemeanscon1 <- data.frame(count=genecountcon1, gene_id=genecount$gene_id) %>% group_by(gene_id) %>% summarise(gm=as.integer(mean(count)))
# genemeanscon2 <- data.frame(count=genecountcon2, gene_id=genecount$gene_id) %>% group_by(gene_id) %>% summarise(gm=as.integer(mean(count)))
# genevarcon1 <- rowVars(allsrr[,1:4])
# genevarcon2 <- rowVars(allsrr[,5:8])

set.seed(1)
## determine isoform switch or differential isoform usage for each gene
geneinfo <- do.call(rbind, lapply(ntrs %>% dplyr::select(ntr) %>% unlist, function(x){
    if (x>1){
        if (runif(1)>noise_level){
            if (runif(1)<0.5){
                return(c(x,T,F,F))
            } else {
                if (runif(1)<0.8){
                    return(c(x,F,T,F))
                } else {
                    return(c(x,F,F,T))
                }
            }
        } else {
            return(c(x,F,F,F))
        }
    } else {
        return(c(x,F,F,F))
    }
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


## simulate fold change  (determine fold change for each gene (not transcript yet))
n<-1
isofc <- c()
fold_change <- do.call(rbind, lapply(ntrs %>% dplyr::select(ntr) %>% unlist, function(x){
    tmp <- geneinfo[n,2:4] 
    type <- which(tmp==1) ### this is indicators of which event is this gene
    
    fc <- sample(seq(2, 5), size=1)

    if (sum(tmp)==0) {
        ### type = 0 : not differentially changing
        isofc <- c(isofc, 1)
        n<<-n+1
        c1 = matrix(rep(1,x))
        c2 = matrix(rep(1,x))
        m = cbind(c1,c2)
        return(m)
    } else if (type==1){
        ### type = 1 : DTE (differentially expression of 1 transcript)
        isofc <- c(isofc, fc)
        c1=matrix(rep(1,x))
        c2=matrix(c(fc,rep(1,x-1)))
        n<<-n+1

        if (runif(1)<0.5){
            m=cbind(c1,c2)
        } else {
            m=cbind(c2,c1)
        }
        
        return(m)
    } else if (type==2) {
        ### type = 2 : DTU == isoform switch (differentially expression of 2 transcript)
        if (x>2){
            m = matrix(c(fc, rep(1, x - 1)))
            c2 = matrix(c(1, fc, rep(1,x-2)))
        } else {
            m = matrix(c(fc,1))
            c2 = matrix(c(1, fc))
        }
        n<<-n+1
        
        # for (i in 2:2) {
        # c2 = sample(m[,1])
        # while (identical(m[,1], c2)) {
        #     c2 <- sample(m[,1])
        # }
        if (runif(1)<0.5){
            m = cbind(m, c2)
        } else {
            m = cbind(c2,m)
        }
        
        isofc <- c(isofc,fc)
        return(m)
    } else if (type==3){
        ### type = 3 : DTEDTU (differentially expression of more than 2 transcripts) 
        m = matrix(c(fc, rep(1, x - 2)))
        for (i in 2:2) {
        c2 = sample(m[,1])
        while (identical(m[,1], c2)) {
            c2 <- sample(m[,1])
        }
        m = cbind(m, c2)
        isofc <- c(isofc,fc)
        n<<-n+1

        if (runif(1)<0.5){
            mm = rbind(m,c(1,fc))
        } else {
            mm = rbind(c(1,fc),m)
        }
        
        return(mm)
    }}
}
))


### simulate isoform ratio using dirichlet distribution

switch_iso <- c()
lastgene <- 1
last_txid <- 0
matrix_of_sums <- lapply(1:nrow(genemeans), function(x){
    ntr <- ntrs$ntr[x]
    genesum <- genemeans$gm[x]
    tmp <- geneinfo[x,2:4]
    type <- which(tmp==1)
    dispgene <- dispd[(last_txid+1):(last_txid+ntr),]$disp
    zero_disp <- which(dispgene==0)

    # switch isoform ratio if the event is isoform switch
    if (ntr>2){
        # if number of transcript > 2
        si <- sample(seq(2, ntr), size=1)
        while (si > ntr){
            si <- sample(seq(2, ntr), size=1)
        }
        alphas <-  c(rep(10, si), rep(1, ntr-si))
        alphas[zero_disp] <- 0.00000001
        originalgeneratio <- rdirichlet(1, alphas)
        if (si==ntr){
            flip_ori_generatio <- sample(originalgeneratio)
        } else {
          flip_ori_generatio <- c(sample(originalgeneratio[1:si]), originalgeneratio[(si+1):ntr])
        }
    } else {
        si <- ntr
        alphas <- rep(10, ntr)
        alphas[zero_disp] <- 0.00000001
        originalgeneratio <- rdirichlet(1, alphas)
        flip_ori_generatio <- 1-originalgeneratio
    }
    
    switch_iso <- c(switch_iso, si)
    flip_ori_generatio[flip_ori_generatio==1] <- 0

    lastgene <<- lastgene+ntr
    if (sum(tmp)==0) {
        c1<-matrix(originalgeneratio)
        c2<-matrix(originalgeneratio)

        m<-cbind(c1,c2)
        return(m)
    } else if (type==1){
        
        c1<-matrix(originalgeneratio)
        c2<-matrix(originalgeneratio)

        m=cbind(c1,c2)
        
        return(m)
    } else if (type==2 | type==3) {
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


### fold change indicating matrix with replicates
fold_change_rep <- do.call(cbind, lapply(1:length(group1n), function(x){
    matrix(rep(as.numeric(t(fold_change[,x])), each = group1n[x]), nrow = nrow(fold_change), byrow = TRUE)
}))


### isoform ratio indicating matrix with replicates
gene_ratio_rep <- do.call(cbind, lapply(1:length(group1n), function(x){
    matrix(rep(as.numeric(t(gene_ratio[,x])), each = group1n[x]), nrow = nrow(gene_ratio), byrow = TRUE)
}))


## expected gene count (gene mean )
baseexp <- matrix(expectedgenecountcon1)
baseexp <- cbind(baseexp, replicate(sum(group1n)-1,baseexp[,1]))
baseexp <- as.data.frame(baseexp)
# baseexp$gene_id <- srr$gene_id
# baseexp <- baseexp %>% group_by(gene_id) %>% mutate(across(colnames(allsrr), sum)) %>% ungroup %>% dplyr::select(-gene_id)

### generate simulated isoform expressions that will be used in negative binomial dis.
final <- baseexp * gene_ratio_rep * fold_change_rep
### calculate dispersion
row.names(allsrr) <- srr$transcript_id
colnames(allsrr) <- lapply(resultfiles, function(x){strsplit(x, "/")[[1]][[1]]}) %>% unlist

### generate simulated isoform expressions for each condition using dispersion from DESeq2 and real gene count
simdf <- do.call(cbind, lapply(1:ncol(final), function(xx){
    tmp1 <- final[,xx]
    
    ## loop for each transcript and apply negative binomial to generate replicates
    eachcol <- do.call(c, lapply(1:length(tmp1), function(xxx){ 
        if (dispd$disp[xxx]==0){
            return(0)
        } else {
            
            size<-1/(dispd$disp[xxx])
            
            if (is.na(size)){
                print(size)
                print(tmp1[xxx])
                print(dispd$disp[xxx]) ####todo check 0 dispersion
            }
            rnbinom(n=1, size=size, mu=tmp1[xxx]) ## 
        }
    }))
    eachcol
}))
simdf[is.na(simdf)] <- 0

#### meta data and groundtruth
colnames(geneinfo) <- c("ntr", "DTE","IS", "DTU")
geneinfo <- as.data.frame(geneinfo)

isoinfo <- genecount %>% 
                dplyr::select(transcript_id, gene_id, iso_ratio) %>% cbind(as.data.frame(mcols(dds))['dispGeneEst'])
isoinfo$isoratioC1 <- gene_ratio[,1]
isoinfo$isoratioC2 <- gene_ratio[,2]
geneinfo$feature_id <- unique(isoinfo$gene_id) %>% unlist
isoinfo$fc <- pmax(fold_change[,1], fold_change[,2])
isoinfo$rowsum <- rowSums(simdf)
isoinfo$DE.status.bin <- ifelse(isoinfo$rowsum==0, 0, ifelse(fold_change[,1] + fold_change[,2]>2, 1, 0))

isoinfo$diffiso <- round(abs(gene_ratio[,2]-gene_ratio[,1]),1)

isoinfo <- isoinfo %>% group_by(gene_id) %>% mutate(new_status=sum(DE.status.bin), n_tx=n()) %>% ungroup
isoinfo <- isoinfo %>% 
                dplyr::mutate(gene_group = ifelse(n_tx>9, ">9", ifelse(n_tx>=5, "5-9", ifelse(n_tx>=2, "2-4", "1"))))


write.table(simdf, paste0(outdir, "/results/simdf.csv", sep="\t"))
write.table(isoinfo, paste0(outdir, sprintf("/results/isoinfo.csv")), sep="\t", row.names = F)
write.table(geneinfo, paste0(outdir, sprintf("/results/geneinfo.csv")), sep="\t", row.names = F)


if (group1n[1]==8){
    sample_name <- c(metadata$sample_id, paste0(metadata$sample_id,"T"))
} else {
    sample_name <- c(metadata$sample_id)
}

#### generate RSEM output for each column using the sample name from the real dataset (or cell ID in single cells)
lapply(1:sum(group1n), function(x){
    tmp <- data.frame(transcript_id=srr$transcript_id, gene_id=srr$gene_id, effective_length=srr$effective_length)
    row.names(tmp) <- tmp$transcript_id
    
    tmp$expected_count <- simdf[,x]
    scale_factors<-simdf[,x] / tmp$effective_length
    scale_factors[is.na(scale_factors) | scale_factors==Inf] <- 0
    tmp_tpm <- (simdf[,x]/tmp$effective_length)*10^6/sum(scale_factors)
    tmp_tpm[is.na(tmp_tpm) | tmp_tpm==Inf] <- 0
    tmp_fpkm <- (simdf[,x]/tmp$effective_length)*10^9/sum(simdf[,x])
    tmp_fpkm[is.na(tmp_fpkm) | tmp_fpkm == Inf]<-0
    tmp$TPM <- tmp_tpm

    tmp$FPKM <- tmp_fpkm 
    sum <- tmp %>% group_by(gene_id) %>% mutate(sumgene=sum(expected_count))
    isopct <- round(tmp$expected_count/sum$sumgene*100, 2)
    isopct[is.na(isopct)]<-0
    tmp$IsoPct <- isopct
    tmp <- tmp[originorder,]
    write.table(tmp, paste0(outdir,"/rsem_sim/", paste0(sample_name,".isoforms.results")[x]), sep="\t", quote=F, row.names=F)
    print("done")
    print(sample_name[x])
})

if (group1n[1]==8){
    newsamples <- c(metadata$sample_id, paste0(metadata$sample_id,"T"))
    write.table(data.frame(sample_id=newsamples, group=c(rep("N", group1n[1]),rep("T",group1n[2]))), paste0(outdir,"/meta.txt"), quote=F, sep="\t",row.names=F)
} else {
    write.table(data.frame(sample_id=metadata$sample_id[1:sum(group1n)], group=c(rep("N", group1n[1]),rep("T",group1n[2]))), paste0(outdir,"/meta.txt"), quote=F, sep="\t",row.names=F)
}

sessionInfo()
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS

# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     

# other attached packages:
#  [1] DESeq2_1.34.0               SummarizedExperiment_1.22.0
#  [3] Biobase_2.52.0              MatrixGenerics_1.4.3       
#  [5] matrixStats_0.61.0          GenomicRanges_1.46.1       
#  [7] GenomeInfoDb_1.30.0         IRanges_2.26.0             
#  [9] S4Vectors_0.30.2            BiocGenerics_0.38.0        
# [11] forcats_0.5.1               stringr_1.4.0              
# [13] dplyr_1.0.6                 purrr_0.3.4                
# [15] readr_1.4.0                 tidyr_1.1.3                
# [17] tibble_3.1.6                ggplot2_3.3.5              
# [19] tidyverse_1.3.1            

# loaded via a namespace (and not attached):
#  [1] bitops_1.0-7           fs_1.5.0               lubridate_1.7.10      
#  [4] bit64_4.0.5            RColorBrewer_1.1-2     httr_1.4.2            
#  [7] tools_4.1.0            backports_1.2.1        utf8_1.2.2            
# [10] R6_2.5.1               DBI_1.1.1              colorspace_2.0-3      
# [13] withr_2.5.0            tidyselect_1.1.2       bit_4.0.4             
# [16] compiler_4.1.0         cli_3.1.0              rvest_1.0.2           
# [19] xml2_1.3.2             DelayedArray_0.18.0    scales_1.1.1          
# [22] genefilter_1.74.1      XVector_0.32.0         pkgconfig_2.0.3       
# [25] dbplyr_2.1.1           fastmap_1.1.0          rlang_1.0.2           
# [28] readxl_1.3.1           rstudioapi_0.13        RSQLite_2.2.11        
# [31] generics_0.1.2         jsonlite_1.7.2         BiocParallel_1.28.3   
# [34] RCurl_1.98-1.6         magrittr_2.0.1         GenomeInfoDbData_1.2.6
# [37] Matrix_1.3-4           Rcpp_1.0.8.2           munsell_0.5.0         
# [40] fansi_1.0.2            lifecycle_1.0.1        stringi_1.7.6         
# [43] zlibbioc_1.38.0        grid_4.1.0             blob_1.2.1            
# [46] crayon_1.5.0           lattice_0.20-44        Biostrings_2.60.2     
# [49] haven_2.4.1            splines_4.1.0          annotate_1.70.0       
# [52] hms_1.0.0              KEGGREST_1.32.0        locfit_1.5-9.5        
# [55] pillar_1.7.0           geneplotter_1.70.0     reprex_2.0.1          
# [58] XML_3.99-0.9           glue_1.6.2             modelr_0.1.8          
# [61] vctrs_0.3.8            png_0.1-7              cellranger_1.1.0      
# [64] gtable_0.3.0           assertthat_0.2.1       cachem_1.0.6          
# [67] xtable_1.8-4           broom_0.7.12           survival_3.2-11       
# [70] AnnotationDbi_1.54.1   memoise_2.0.1          ellipsis_0.3.2     