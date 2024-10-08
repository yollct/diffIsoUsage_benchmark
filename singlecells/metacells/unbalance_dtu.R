library(DTUrtle)
library(satuRn)
library(tidyverse)
library(biomaRt)
library(satuRn)
library(SummarizedExperiment)
library(tximport)
library(DEXSeq)
library(DRIMSeq)

library(edgeR)
source('/nfs/proj/is_benchmark/singlecells/testDTU.R')
biocpar <- BiocParallel::MulticoreParam(4)
outdir <- "/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/"

reps <- c(10, 50, 100, 200, 500)#n of meta cells
allfolders <- list.files(outdir, pattern="sim_data_unbalance")
allfolders <- allfolders[grepl("0$",allfolders)] ### only noise = 0


human <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=107)

resgene <- do.call(c, lapply(allfolders, function(thisfolder){
    new_count <- readRDS(sprintf("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/%s/sim_unbalance_2ct_rep%s.rds", thisfolder, 7))

    # get effective length
    ensembl_list <- row.names(new_count)
    gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_transcript_id", "transcript_start","transcript_end"), filters="ensembl_transcript_id", values=ensembl_list, mart=human)
    gene_coords$size=gene_coords$transcript_end - gene_coords$transcript_start
    gene_coords <- gene_coords %>% dplyr::select(ensembl_transcript_id, size)

    gene_length <- data.frame(transcript_id = row.names(new_count))
    gene_length <- left_join(gene_length, gene_coords, by=c("transcript_id"="ensembl_transcript_id")) 
    gene_length <- gene_length %>% group_by(transcript_id) %>% summarise(size=max(size)) %>% as.data.frame()
    sum(is.na(gene_length$size))

    # Function to write each column of the count matrix to a Salmon-like quant.sf file
    write_salmon_quant_files <- function(count_matrix) {
    for (sample in colnames(count_matrix)) {
        # Extract counts for the current sample
        counts <- count_matrix[, sample]
        
        # Create a data frame similar to Salmon's quant.sf
        quant_sf <- tibble(
        Name = rownames(count_matrix),
        Length = NA,  # Assuming Length is not available
        EffectiveLength = gene_length$size,  # Assuming EffectiveLength is not available
        TPM = (counts / gene_length$size) / sum(counts / gene_length$size) * 1e6,  # Assuming TPM is not calculated here
        NumReads = counts
        )
        dir.create(paste0(outdir, "unbalance", 7))
        dir.create(paste0(outdir, "unbalance",7,"/",sample))
        # File path for the output
        file_path <- paste0(outdir, "unbalance",7,"/", sample, "/quant.sf")
        
        # Write the data frame to a file
        write_tsv(quant_sf, file_path, col_names = TRUE)
    }
    }

    # Run the function with your count matrix
    write_salmon_quant_files(new_count)
    files <- Sys.glob(sprintf("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/unbalance%s/*/quant.sf", 7))
    names(files) <- gsub(".*/","",gsub("/quant.sf","",files))
    cnt <- import_counts(files = files, type = "salmon")
    metadata <- read.csv(sprintf("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/%s/sim_unbalance_2ct_rep%s_metadata.csv", thisfolder, 7), sep="\t")
    metadata$group <- paste0("celltype", metadata$group)
    # get effective length
    meta <- data.frame("sample_id"=colnames(cnt), "group"=factor(metadata$group, levels=unique(metadata$group)), 
                    stringsAsFactors = FALSE)

    res_g <- lapply(reps, function(r){
        print(r)

        tx2gene <- import_gtf(gtf_file = "/nfs/proj/is_benchmark/singlecells/Homo_sapiens.GRCh38.107.gtf.gz")
        tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_id", "gene_id"))

        

        # cnt <- import_counts(files = files, type = "salmon")

        # txinfo <- genetx %>% dplyr::select(gene_id, transcript_id) %>% dplyr::rename(isoform_id = "transcript_id")
        # row.names(txinfo) <- txinfo$isoform_id
        # txinfo <- txinfo[which(txinfo$isoform_id %in% row.names(cts$counts)),]
        # txinfo <- subset(txinfo, duplicated(txinfo$gene_id) | duplicated(txinfo$gene_id, fromLast = TRUE))

        # cnt <- cts$counts[which(row.names(cts$counts) %in% txinfo$isoform_id),]
        #txinfo1 <- txinfo[which(row.names(cts$counts) %in% txinfo$isoform_id),]
        # meta$sample_name <- meta$sample_id
        # row.names(meta) <- NULL

        # Modified code
        subcnt <- do.call(cbind, lapply(seq(1, 3), function(x) {
            # Initialize an empty list to store sub-counts for each cell type
            subcnt_list <- list()
            
            # Loop over each cell type
            for (ct in 1:7) {
                # Sample 'r' samples from each cell type
                samples <- sample(meta[meta$group == sprintf("celltype%d", ct), ]$sample_id, r, replace=T)
                
                # Compute the row sums of counts for the selected samples
                ct_counts <- rowSums(cnt[, samples])
                
                # Add the computed counts to the list and name it
                subcnt_list[[ct]] <- ct_counts
                names(subcnt_list)[ct] <- sprintf("metacell%d_%s", ct, x)
            }
            
            # Combine counts from all cell types for this iteration
            do.call(cbind, subcnt_list)
        }))

        pd <- data.frame("sample_id"=colnames(subcnt))
        pd$group <- lapply(pd$sample_id, function(x){strsplit(x, "_")[[1]][1]}) %>% unlist
        pd$rep <- lapply(pd$sample_id, function(x){strsplit(x, "_")[[1]][2]}) %>% unlist

        allcomb <- apply(combn(unique(pd$group), 2),2,paste0,collapse='_')

        dturtle_g <- do.call(rbind, lapply(allcomb, function(x){
            print(x)
            cell1 <- strsplit(x, "_")[[1]][1]
            cell2 <- strsplit(x, "_")[[1]][2]            
            dturtle <- run_drimseq(counts = subcnt, tx2gene = tx2gene, pd=pd, id_col = "sample_id",
                                cond_col = "group", cond_levels = c(cell1, cell2), filtering_strategy = "bulk", 
                                BPPARAM = biocpar)
            dturtle_1 <- posthoc_and_stager(dturtle = dturtle, ofdr = 1, posthoc = 0)
            dturtle_g <- data.frame(feature_id=dturtle_1$FDR_table$geneID, fdr=dturtle_1$FDR_table$gene, tool="DTUrtle", nmetacells=r, comb=x) %>% unique
            dturtle_tx <- data.frame(feature_id=dturtle_1$FDR_table$txID, fdr=dturtle_1$FDR_table$transcript, tool="DTUrtle", nmetacells=r, comb=x)
            dturtle_g
        }))

        dturtle_g <- dturtle_g %>% reframe(fdr = min(fdr), tool=tool, nmetacells=r, .by="feature_id") %>% unique
        
        ## saturn
        saturn_g <- do.call(rbind, lapply(allcomb, function(x){
            print(x)
            cell1 <- strsplit(x, "_")[[1]][1]
            cell2 <- strsplit(x, "_")[[1]][2]   

            ## saturn
            filter_edgeR <- filterByExpr(subcnt,
                                design = NULL,
                                group = pd$group,
                                lib.size = NULL,
                                min.count = 1,
                                min.prop = 0.1)                    

            thiscts <- subcnt[filter_edgeR,]
            tx2gene <- tx2gene %>% dplyr::rename("isoform_id"="transcript_id")
            tx2gene <-  tx2gene[which(tx2gene$isoform_id %in% row.names(thiscts)),]
            tx2gene <- subset(tx2gene, 
                        duplicated(tx2gene$gene_name) | duplicated(tx2gene$gene_name, fromLast = TRUE))
            thiscts <- thiscts[which(row.names(thiscts) %in% tx2gene$isoform_id), ]
            tx2gene <- tx2gene[match(row.names(thiscts), tx2gene$isoform_id),]
            

            pd$sample_name <- pd$sample_id
            row.names(pd) <- NULL
            sumExp <- SummarizedExperiment::SummarizedExperiment(
                assays = list(counts = thiscts),
                colData = pd,
                rowData = tx2gene
            )
            metadata(sumExp)$formula <- ~ 0 + as.factor(colData(sumExp)$group)
            print("fitting DTU")
            sumExp <- satuRn::fitDTU(
                        object = sumExp,
                        formula = ~ 0 + group,
                        parallel = FALSE,
                        BPPARAM = BiocParallel::bpparam(),
                        verbose = TRUE
                    )
            
            group <- as.factor(pd$group)
            design <- model.matrix(~0+group)
            colnames(design) <- levels(group)
            y <- paste0(cell1,"-",cell2)
            print(y)
            L <- limma::makeContrasts(contrasts=y, levels =c(cell1, cell2))

            print(sumExp)
            print("testing DTU")
            sumExp <- testDTU(
                        object = sumExp,
                        contrasts = L,
                        diagplot1 = FALSE,
                        diagplot2 = FALSE,
                        sort = FALSE
                    )

            print("joining results")
            res <- rowData(sumExp)[[paste0("fitDTUResult_", y)]]
            res$isoform_id <- row.names(res)
            res <- full_join(res, tx2gene, by=c("isoform_id"))
            res$isoform_id <- lapply(res$isoform_id, function(x){strsplit(x, "[.]")[[1]][1]}) %>% unlist()
            
            saturn_g <- res %>% dplyr::select(gene_id, regular_FDR) %>% group_by(gene_id) %>% summarise(saturn=min(regular_FDR)) %>% dplyr::rename("feature_id"="gene_id", "fdr"="saturn") 
            saturn_g$tool <- "satuRn"
            saturn_g$nmetacells <- r
            saturn_g$comb <- x
            
            saturn_g <- data.frame(feature_id=res$gene_id, fdr=res$regular_FDR, tool="satuRn", nmetacells=r, comb=x)
        }))

        saturn_g <- saturn_g %>% reframe(fdr = min(fdr), tool=tool, nmetacells=r, .by="feature_id") %>% unique
        print("done")

        ##limmaDS
        ## limmads 
        limma_g <- do.call(rbind, lapply(allcomb, function(x){
            print(x)
            print("limma")
            cell1 <- strsplit(x, "_")[[1]][1]
            cell2 <- strsplit(x, "_")[[1]][2]   

            thispd <- pd %>% dplyr::filter(group %in% c(cell1, cell2))
            thiscnt <- subcnt[,thispd$sample_id]

            thistx <- tx2gene[tx2gene$transcript_id %in% row.names(thiscnt),]

            dge <- DGEList(counts=thiscnt,genes=thistx$transcript_id)
            dge$genes$GeneID <- thistx$gene_id
            A <- rowSums(dge$counts)
            dge <- dge[A>1,, keep.lib.sizes=F]
            dge <- calcNormFactors(dge)
            design <- model.matrix(~ group, data = thispd)
            v <- voom(dge, design, plot=F)
            fit <- lmFit(v, design)

            ex <- diffSplice(fit, geneid="GeneID")
            res <- topSplice(ex, n=nrow(dge$genes), coef=2, test="t")

            res1 <- res %>% dplyr::select(GeneID, FDR) %>% group_by(GeneID)  %>% summarise(FDR=min(FDR)) %>% na.omit() %>% unique()
            limma_g <- res1 %>% dplyr::rename("feature_id"="GeneID", 'fdr'="FDR")
            limma_g$tool <- "LimmaDS"
            limma_g$nmetacells <- r
            limma_g$comb <- x

            limma_g
        }))
        limma_g <- limma_g %>% reframe(fdr = min(fdr), tool=tool, nmetacells=r, .by="feature_id") %>% unique
        print("done")

        dx_g <- do.call(rbind, lapply(allcomb, function(x){
            print(x)
            print("dexseq")
            cell1 <- strsplit(x, "_")[[1]][1]
            cell2 <- strsplit(x, "_")[[1]][2]   
            ##DEXSeq 
            thistx <- tx2gene[tx2gene$transcript_id %in% row.names(subcnt), ]
            subcnt <- data.frame(subcnt)
            subcnt$gene_id <- thistx$gene_id
            subcnt$feature_id <- thistx$transcript_id
            # subcnt <- subcnt %>% dplyr::mutate(gene_id = thistx$gene_id, feature_id=thistx$transcript_id)

            d <- dmDSdata(counts=subcnt, samples=pd)
            dxd <- DEXSeqDataSet(countData=round(as.matrix(DRIMSeq::counts(d)[,-c(1:2)])), 
                                sampleData=DRIMSeq::samples(d),
                                design=~sample+exon+group:exon,
                                featureID=counts(d)$feature_id,
                                groupID=counts(d)$gene_id)

            dxd <- estimateSizeFactors(dxd)
            dxd <- estimateDispersions(dxd, BPPARAM = biocpar, quiet=T)
            dxd <- testForDEU(dxd, BPPARAM = biocpar, reducedModel = ~sample+exon)
            dxd = estimateExonFoldChanges( dxd, fitExpToVar="group", BPPARAM = biocpar)
            dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
            qval <- perGeneQValue(dxr)
            dxr <- as.data.frame(dxr)
            qval[qval==0] <- min(qval[qval !=0]) * 0.1
            dxr_g <- data.frame(feature_id=names(qval), fdr=qval)
            dxr_g$tool <- "DEXSeq"
            dxr_g$nmetacells <- r
            dxr_g$comb <- x
            dxr_g
        }))
        print(dx_g %>% head)
        dx_g <- dx_g %>% reframe(fdr = min(fdr), tool=tool, nmetacells=r, .by="feature_id") %>% unique
        print("done")
        # unlink(paste0(outdir, "unbalance", 7), recursive = TRUE)
    
        res_g <- do.call("rbind", list(dturtle_g, saturn_g, limma_g, dx_g))
        res_g$folder <- thisfolder

        return(res_g)
        })
        unlink(paste0(outdir, "unbalance", 7), recursive = TRUE)
        res_g
    })
    

 )
# print(resgene %>% head)
resgene <- do.call(rbind, resgene)

tools <- c("DTUrtle", "satuRn", "LimmaDS", "DEXSeq")
testpos <- lapply(allfolders, function(thisfolder){
    lapply(reps, function(r){
    lapply(tools, function(x){
        resgenefdr <- resgene %>% dplyr::filter(folder==thisfolder & tool==x & nmetacells==r & fdr < 0.05)
        resgenefdr$feature_id
    })})
}) 
testpos<- unlist(testpos, recursive=F)
testpos <- unlist(testpos, recursive=F)
names(testpos) <- Map(paste0, base::rep(allfolders, each=length(reps)*length(tools)), tools, base::rep(reps, each=4), sep="")
saveRDS(testpos, "testpos_genes.rds")

#### calculate precision and recall
truepos <- lapply(allfolders, function(thisfolder){
    lapply(reps, FUN=function(r){
    # Assuming count_matrix is your single-cell count matrix
    # Replace this with your actual matri
    print("run")
    print(r)
    new_geneinfo <- read.csv(sprintf("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/%s/sim_unbalance_2ct_rep%s_geneinfo.csv", thisfolder, 7), sep="\t")
    
    truepos <- new_geneinfo[new_geneinfo$DTU==1,]$feature_id 

    # truepos <- new_geneinfo[new_geneinfo$DE.status.bin==1,]$transcript_id
    unique(truepos)
})})
truepos <- unlist(truepos, recursive=F)
names(truepos) <- Map(paste0, base::rep(allfolders, each=length(reps)), reps)

#calculate precision recall
cal_pre_re <- function(testpos, truepos) {
    xx <- 1
    mets <- c("DTUrtle", "satuRn", "LimmaDS", "DEXSeq")

    out2 <- do.call(c, lapply(allfolders, function(thisfolder) {
        out <- do.call(c, lapply(mets, function(method) {
            lapply(reps, function(r) {
                print(paste0(method, r))
                
                testgenes <- testpos[[paste0(thisfolder, method, r)]]
                truegenes <- truepos[[paste0(thisfolder, as.character(r))]]

                print(testgenes %>% head)
                tp <- sum(testgenes %in% truegenes)
                fp <- sum(!testgenes %in% truegenes)
                fn <- sum(!truegenes %in% testgenes)

                pp <- length(unique(testgenes))
                print(tp)
                print(pp)
                precisions <- tp / pp

                p <- length(truegenes)
                recall <- tp / p

                f1 <- 2 * tp / (2 * tp + fp + fn)

                xx <- xx + 1

                outout <- data.frame(c(precisions, recall, method, r, thisfolder))
                t(outout)
            })
        }))
        out
    }))

    out2
}


outres <- cal_pre_re(testpos, truepos)
final <- do.call(rbind, outres)
colnames(final)<-c("precision", "recall", "tool", "nmetacells", "folder")
row.names(final) <- seq(1, nrow(final))
final <- data.frame(final)
final

final_long <- pivot_longer(final, c("precision", "recall"), names_to="measure", values_to="values")
final_long$values <- as.numeric(final_long$values)
final_long$nmetacells <- as.numeric(final_long$nmetacells)
final_long <- final_long %>% group_by(tool, nmetacells, measure) %>% summarise(mean=mean(values), sd=sd(values))

write.table(final_long, "/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/metacell_unbalance_finallong.csv")

final_long <- read.csv("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/metacell_unbalance_finallong.csv", sep=" ")
supp.labs <- c("Precision", "Recall")
names(supp.labs) <- c("precision", "recall")


png("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/metacell_0_unbalance_prere.png", res=300, width=3000, height=1000)
ggplot(final_long, aes(x=nmetacells, y=mean, color=tool))+
    geom_point()+
    geom_line()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color=tool), width=.2,
                 position=position_dodge(0.05))+
    facet_grid(.~measure, labeller=labeller(measure=supp.labs))+
    scale_color_manual(values=c("DTUrtle"='#999999','satuRn'='#E69F00','DEXSeq'='#0077b6','LimmaDS'='#e35d6a'))+
    theme_bw()+xlab("Number of cells per group") +ylab("Precision/Recall")
dev.off()

