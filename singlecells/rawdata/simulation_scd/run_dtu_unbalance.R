
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
allfolders <- allfolders[grepl("0$",allfolders)] ### only noise = 0

human <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=107)


resgene <- do.call(c, lapply(allfolders, function(thisfolder){
    lapply(reps, function(r){
    # Assuming count_matrix is your single-cell count matrix
    # Replace this with your actual matri
    print("run")
    print(r)
    new_count <- readRDS(sprintf("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/%s/sim_unbalance_2ct_rep%s.rds", thisfolder, r))

    metadata <- read.csv(sprintf("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/%s/sim_unbalance_2ct_rep%s_metadata.csv", thisfolder, r), sep="\t")
    metadata$group <- paste0("celltype", metadata$group)
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
        dir.create(paste0(outdir, "unbalance", r))
        dir.create(paste0(outdir, "unbalance",r,"/",sample))
        # File path for the output
        file_path <- paste0(outdir, "unbalance",r,"/", sample, "/quant.sf")
        
        # Write the data frame to a file
        write_tsv(quant_sf, file_path, col_names = TRUE)
    }
    }

    # Run the function with your count matrix
    write_salmon_quant_files(new_count)

    tx2gene <- import_gtf(gtf_file = "/nfs/proj/is_benchmark/singlecells/Homo_sapiens.GRCh38.107.gtf.gz")
    tx2gene <- move_columns_to_front(df = tx2gene, columns = c("transcript_id", "gene_id"))

    files <- Sys.glob(sprintf("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/unbalance%s/*/quant.sf", r))
    names(files) <- gsub(".*/","",gsub("/quant.sf","",files))
    cts <- import_counts(files = files, type = "salmon")

    # cts <- combine_to_matrix(tx_list = cts)

    pd <- data.frame("sample_id"=colnames(cts), "group"=factor(metadata$group, levels=unique(metadata$group)), 
                    stringsAsFactors = FALSE)

    allcomb <- apply(combn(unique(pd$group), 2),2,paste0,collapse='_')

    dturtle_g <- do.call(rbind, lapply(allcomb, function(x){
        cell1 <- unique(pd$group)[as.integer(strsplit(x, "_")[[1]][1])]
        cell2 <- unique(pd$group)[as.integer(strsplit(x, "_")[[1]][2])]
        print(cell1)
        print(cell2)
        dturtle <- run_drimseq(counts = cts, tx2gene = tx2gene, pd=pd, id_col = "sample_id",
                            cond_col = "group", cond_levels = c(cell1, cell2), filtering_strategy = "sc", 
                            BPPARAM = biocpar)
        dturtle_1 <- posthoc_and_stager(dturtle = dturtle, ofdr = 1, posthoc = 0)
        dturtle_g <- data.frame(feature_id=dturtle_1$FDR_table$geneID, fdr=dturtle_1$FDR_table$gene, tool="DTUrtle", rep=r, comb=x) %>% unique
        dturtle_tx <- data.frame(feature_id=dturtle_1$FDR_table$txID, fdr=dturtle_1$FDR_table$transcript, tool="DTUrtle", rep=r, comb=x)
        dturtle_g
    }))

    dturtle_g <- dturtle_g %>% reframe(fdr = min(fdr), tool=tool, rep=r, .by="feature_id") %>% unique
    
    saturn_g <- do.call(rbind, lapply(allcomb, function(x){
        cell1 <- unique(pd$group)[as.integer(strsplit(x, "_")[[1]][1])]
        cell2 <- unique(pd$group)[as.integer(strsplit(x, "_")[[1]][2])]
        print(cell1)
        print(cell2)

        ## saturn
        filter_edgeR <- filterByExpr(cts,
                            design = NULL,
                            group = pd$group,
                            lib.size = NULL,
                            min.count = 1,
                            min.total.count = 5,
                            large.n = 3,
                            min.prop = 0.1)                    

        thiscts <- cts[filter_edgeR,]
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
        saturn_g$rep <- r
        saturn_g$comb <- x
        
        saturn_g <- data.frame(feature_id=res$gene_id, fdr=res$regular_FDR, tool="satuRn", rep=r, comb=x)
    }))

    saturn_g <- saturn_g %>% reframe(fdr = min(fdr), tool=tool, rep=r, .by="feature_id") %>% unique
    print("done")

    unlink(paste0(outdir, "unbalance", r), recursive = TRUE)
    
    res_g <- rbind(dturtle_g, saturn_g)
    res_g$folder <- thisfolder
    # res_tx <- rbind(dturtle_tx, saturn_tx)
    return(res_g)
    })})  
)
 
resgene <- do.call(rbind, resgene)
print(resgene %>% head)
tools <- c("DTUrtle", "satuRn")
testpos <- lapply(allfolders, function(thisfolder){
    lapply(reps, function(r){
    lapply(tools, function(x){
        resgenefdr <- resgene %>% dplyr::filter(folder==thisfolder & tool==x & rep==r & fdr < 0.05)
        resgenefdr$feature_id
    })})
}) 
testpos<- unlist(testpos, recursive=F)
testpos <- unlist(testpos, recursive=F)
names(testpos) <- Map(paste0, base::rep(allfolders, each=12), tools, base::rep(reps, each=2), sep="")
saveRDS(testpos, "testpos_genes.rds")
#### calculate precision and recall
truepos <- lapply(allfolders, function(thisfolder){
    lapply(reps, FUN=function(r){
    # Assuming count_matrix is your single-cell count matrix
    # Replace this with your actual matri
    print("run")
    print(r)
    new_geneinfo <- read.csv(sprintf("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/%s/sim_unbalance_2ct_rep%s_geneinfo.csv", thisfolder, r), sep="\t")
    
    truepos <- new_geneinfo[new_geneinfo$DTU==1,]$feature_id 
    # truepos <- new_geneinfo[new_geneinfo$DE.status.bin==1,]$transcript_id
    unique(truepos)
})})
truepos <- unlist(truepos, recursive=F)
names(truepos) <- Map(paste0, base::rep(allfolders, each=6), reps)
saveRDS(truepos, "truepos_genes.rds")

#calculate precision recall
#calculate precision recall
cal_pre_re <- function(testpos, truepos) {
    xx <- 1
    mets <- c("DTUrtle", "satuRn")

    out2 <- do.call(c, lapply(allfolders, function(thisfolder) {
        out <- do.call(c, lapply(mets, function(method) {
            lapply(reps, function(r) {
                print(paste0(method, r))
                
                testgenes <- testpos[[paste0(thisfolder, method, r)]]
                truegenes <- truepos[[paste0(thisfolder, as.character(r))]]

                tp <- sum(testgenes %in% truegenes)
                fp <- sum(!testgenes %in% truegenes)
                fn <- sum(!truegenes %in% testgenes)

                pp <- length(unique(testgenes))
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
colnames(final)<-c("precision", "recall", "tool", "reps", "folder")
row.names(final) <- seq(1, nrow(final))
final <- data.frame(final)
final
write.table(final, "singlecell_0_unbalance_final.csv", row.names=F, quote=F)

final_long <- pivot_longer(final, c("precision", "recall"), names_to="measure", values_to="values")
final_long$values <- as.numeric(final_long$values)
final_long$reps <- as.numeric(final_long$reps)
final_long <- final_long %>% group_by(tool, reps, measure) %>% summarise(mean=mean(values), sd=sd(values))

write.table(final_long, "singlecell_unbalance_finallong.csv", row.names=F, quote=F)
final_long <- read.csv("/nfs/proj/is_benchmark/singlecells/rawdata/simulation_scd/singlecell_unbalance_finallong.csv", sep=" ")
supp.labs <- c("Precision", "Recall")
names(supp.labs) <- c("precision", "recall")
png("singlecell_0_unbalance_prere.png", res=300, width=3000, height=1000)
ggplot(final_long, aes(x=reps, y=mean, color=tool))+
    geom_point()+
    geom_line()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                 position=position_dodge(0.05))+
    facet_grid(.~measure, labeller = labeller(measure=supp.labs))+
    scale_color_manual(values=c("DTUrtle"='#999999','satuRn'='#E69F00'))+
    theme_bw()+xlab("Number of cell types")+ylab("Precision/Recall")+
    theme(text = element_text(size = 15),
        axis.text.x = element_text(size=15, hjust = 1),
        axis.text.y = element_text(size=15, hjust = 1)) 
dev.off()

mets <- c("DTUrtle", "satuRn")
list_false_negative <- lapply(reps, function(r){
        this <- lapply(mets, function(method){
            testgenes <- testpos[[paste0(method,r)]]
            truegenes <- truepos[[as.character(r)]]

            truegenes[!truegenes %in% testgenes]
    })
    names(this) <- mets
    this
})
names(list_false_negative) <- reps
saveRDS(list_false_negative, 'listfalseneg.rds')

list_false_positive <- lapply(reps, function(r){
        this <- lapply(mets, function(method){
            testgenes <- testpos[[paste0(method,r)]]
            truegenes <- truepos[[as.character(r)]]

            testgenes[!testgenes %in% truegenes]
    })
    names(this) <- mets
    this
})
names(list_false_positive) <- reps
saveRDS(list_false_negative, 'listfalsepos.rds')
