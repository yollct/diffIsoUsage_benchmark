library(tidyverse)
library(ggplot2)
library(UpSetR)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(GGally)
library(ggradar)
library(fmsb)
library(gridExtra)
library(scales)
source("/nfs/proj/is_benchmark/plot_functions/new_functions.R")

path <- "/nfs/proj/is_benchmark"
outdir <- "/nfs/scratch/chit/new_simulations/%s_50_%s_%s%s"
noises <- c("_0", "_0.1","_0.5", "_0.7", "_0.9")

thisrep <- c('4','8')
seqtype <- 'pair'
whichrep <- 'r1'

thresholds <- data.frame(tools=c("iso_ktsp","drimseq","dexseq","dturtle","seqGSEA","cuffdiff","junctionseq","saturn","drimseq_stageR", "dexseq_stageR","saturn_stageR", "DSGseq", "nbsplice", "LimmaDS", "edgeR", "BANDITS"),
                        thres=c(0.5,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05, 5, 0.05, 0.05, 0.05, 0.05))
row.names(thresholds) <- thresholds$tools

get_overall_tpfp <- function(split){
    overalldf <- do.call(rbind, lapply(whichrep, function(R){
            do.call(rbind, lapply(seqtype, function(seq){
            do.call(rbind, lapply(thisrep, function(rep){
                do.call(rbind, lapply(noises, function(noise){
                    print(paste0(sprintf(outdir, seq, rep, R, noise)))
                    dirr <- sprintf(outdir, seq, rep, R, noise)
                    
                    truthfiles_g <- read.csv(paste0(dirr, "/results/truthtable_gene.csv"), sep="\t") ### change to result
                    truthfiles_tx <- read.csv(paste0(dirr, "/results/truthtable_tx.csv"), sep="\t")
                    
                    if (is.null(split)){
                        
                        # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                    } else {
                        if (split!="events"){
                        print("not consider DTE")
                        # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                        }
                    }

                    diffiso <- read.csv(paste0(dirr, "/results/isoinfo.csv"), sep="\t")
                    isorat <- diffiso %>% group_by(gene_id) %>% summarise(diffiso=max(diffiso))
                    truthfiles_g$diffiso <- isorat$diffiso
                    truthfiles_g <- truthfiles_g %>% dplyr::mutate(diffiso=ifelse(diffiso>0.8, ">0.8", ifelse(diffiso>=0.6, "0.6-0.8", ifelse(diffiso>=0.5, "0.5-0.6", ifelse(diffiso>=0.4, "0.4-0.5", ifelse(diffiso>=0.3, "0.3-0.4", ifelse(diffiso>=0.2, "0.2-0.3", ifelse(diffiso>=0.1, "0.1-0.2", "0-0.1")))))))) 
                    print(truthfiles_g %>% nrow)
                    # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    # truthfiles_tx <- truthfiles_tx %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    
                    do.call(rbind, lapply(c("kal", "salmon", "rsem"), function(quant){
                        df_g <- read.csv(paste0(dirr, "/results/", sprintf("%s_res_tx_N_T.txt", quant)), sep="\t")
                        df_g <- df_g %>% dplyr::select(-any_of("BayesDRIMSeq"))
                        thisthres <- thresholds[colnames(df_g)[2:ncol(df_g)],]$thres
                        df_g <- df_g %>% group_by(feature_id) %>% summarise(across(colnames(select(df_g, -feature_id)), .fns = min), .groups = "keep")

                        df_g <- df_g[grepl("ENST", df_g$feature_id),]

                        this_pr <- cat_tp_fp(df_g, truthfiles_tx)                        
                        
                        this_pr$quant_tool <- quant
                        this_pr$rep <- paste(as.character(rep), 'replicates')
                        this_pr$seq <- seq
                        this_pr$sim <- R
                        this_pr$con <- paste0(quant, "_", paste(as.character(rep), 'replicates'))
                        this_pr$noise <- gsub("_", "", noise)
                        rm(df_g)
                        this_pr
                    }))
                }))
            }))
        }))
    }))
    overalldf$groupcount <- lapply(overalldf$groupcount, function(x){as.numeric(x)}) %>% unlist
    overalldf
}

cat_tp_fp <- function(df_g, truthfiles_tx){
    # Define the method names (all columns except the first)
    methods <- colnames(df_g)[-1]

    # Retrieve threshold values for each method in the same order.
    # This assumes thresholds is indexed by method name and has a column "thres".
    threshold_values <- thresholds[methods, "thres"]

    # Use sapply to iterate over each method index,
    # process the column, and return a logical vector.
    out <- sapply(seq_along(methods), function(i) {
      method <- methods[i]

      t_val  <- threshold_values[i]
      
      # Get the vector for this method
      col_vals <- df_g[[method]]
      print(method)
      # Replace NA values with 1 if threshold equals 0.05, else with 0.
      # This is vectorized.
      col_vals[is.na(col_vals)] <- if (t_val == 0.05) 1 else 0
      
      # Return TRUE/FALSE based on the condition:
      # - if t_val equals 0.05, check if col_vals < t_val;
      # - otherwise, check if col_vals > t_val.
      if (t_val == 0.05) {
        col_vals < t_val
      } else {
        col_vals > t_val
      }
    })
    colnames(out) <- colnames(df_g[,-1])
    out <- as.data.frame(out)
    out$feature_id <- df_g$feature_id

    joined_res <- full_join(truthfiles_tx %>% dplyr::select(feature_id, gene_id, status, events), out, by=c("feature_id"))
    print(sum(joined_res$iso_ktsp))

    joined_long <- joined_res %>%
      pivot_longer(cols = all_of(methods), 
                  names_to = "tool", 
                  values_to = "detected")
    
    joined_long$detected[is.na(joined_long$detected)] <- FALSE

    gene_pp <- joined_res %>% dplyr::select(gene_id, status) %>% group_by(gene_id) %>% summarise(PP=sum(status))
    tpr <- joined_res %>% group_by(gene_id) %>% group_by(events, status) %>% summarise(P=sum(status)) %>% filter(status==1)

    joined_long$status <- as.numeric(joined_long$status)

    joined_tp <- joined_long %>%
      group_by(gene_id, tool, events) %>%
      summarise(
        TP = sum(status == 1 & detected),
        FP = sum(status == 0 & detected),
        .groups = "drop"
      )
    
    joined_tp <- left_join(joined_tp, gene_pp, by=c("gene_id"))

    joined_tp <- joined_tp %>% filter(PP!=0)
    

    joined_tp <- joined_tp %>%
      mutate(category = case_when(
        TP == PP & TP!=0 & FP == 0 ~ "All detected",
        TP > 0 & FP == 0 ~ "Only TP",
        TP > 0 & FP > 0  ~ "TP and FP",
        TP == 0 & FP > 0 ~ "Only FP",
        TRUE             ~ "None"   # Covers the case where TP == 0 and FP == 0
      ))
    
    
    joined_tp <- joined_tp %>% group_by(tool, events, category) %>% summarise(groupcount=n())
    joined_tp <- left_join(joined_tp, tpr, by=c("events"))

    return(joined_tp)
}
rename_if_exists <- function(df, rename_map) {
  for (old_name in names(rename_map)) {
    if (old_name %in% names(df)) {
      names(df)[names(df) == old_name] <- rename_map[[old_name]]
    }
  }
  return(df)
}

newnamelist <- list("dexseq"="DEXSeq", "drimseq"="DRIMSeq", "iso_ktsp"="isoKTSP", "saturn"="satuRn", "dturtle"="DTUrtle", "nbsplice"="NBSplice")

tpfp_overalldf <- get_overall_tpfp(split=NULL)

tpfp_overalldf$normgroupcount <- ifelse(tpfp_overalldf$groupcount==0, 0, tpfp_overalldf$groupcount/tpfp_overalldf$P)
tpfp_overalldf <- tpfp_overalldf %>% filter(category!="None" & status!=0)

plot_list <- c()
plot_name <- c()
n<-1
for (spl in sort(unique(tpfp_overalldf$category))){
    for (rps in c("DTE", "DTU", "IS")){
        data <- tpfp_overalldf %>% filter(quant_tool=="salmon"&rep=="4 replicates"&category==spl&events==rps) %>% dplyr::select(events, noise, normgroupcount, tool) %>% pivot_wider(names_from = tool, values_from = normgroupcount, id_cols = noise) %>% as.data.frame()
        
        if (nrow(data)==0){
          
          cols <- colnames(data)
          data <- as.data.frame(matrix(0, nrow = 5, ncol = 4))
          colnames(data) <- c("noise", "BANDITS", "nbsplice", "satuRn")

        }
        data <- rbind(rep(0.2,ncol(data)) , rep(0,ncol(data)) , data)
        data <- data %>% select(-noise) 
        data[is.na(data)] <- 0
        
        data <- rename_if_exists(data, newnamelist)

        plot_list[[n]] <- data
        plot_name[[n]] <- paste0(rps,"_",spl)
        n<-n+1
 }
}


colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9), "#008000" , "black" )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4),"#008000" , "black" )

png(sprintf("./tx_noise_fig/tpfp_cat_radar.png"), width=5000, height=4500, res=300)
nrows<-4
ncols<-3
par(mfrow = c(nrows, ncols), mar = c(1, 2, 1, 2), oma = c(4, 6, 4, 6)) # 4 rows, 2 columns

for (i in seq_along(plot_list)) {
radarchart(plot_list[[i]], axistype=1 , 
        #custom polygon
        pcol=colors_border ,  plwd=4 , plty=1,
        #custom the grid
        cglcol="grey", cglty=1, axislabcol="#5E5E5E", caxislabels=seq(0,0.2,0.05), cglwd=0.8,
        #custom labels
        vlcex=2, calcex=1.5
        
)  
}
row_labels <- sort(unique(tpfp_overalldf$category)) # Adjust as needed
for (i in 1:nrows) {
  mtext(row_labels[i], side = 2, at = (nrows - i + 0.5)/nrows, srt=90, outer = TRUE, line = 1, cex=1.8)
}

# Add column labels
col_labels <- c("S1", "S2", "S3")   # Adjust as needed
for (i in 1:ncols) {
  mtext(col_labels[i], side = 3, at = (i - 0.5)/ncols, las = 1, outer = TRUE, line = 1, cex=1.8)
}

dev.off()

df_long_counts <- joined_tp %>%
  pivot_longer(
    cols = c(TP, FP),
    names_to = "Metric",
    values_to = "Count"
  )

library(ggplot2)
png(sprintf("./tx_noise_fig/tx_gene_tpfp.png"), width=3500, height=2000, res=300)
ggplot(df_long_counts %>% filter(Count!=0), aes(x = Count, fill = Metric)) +
  geom_bar(position = "dodge") +
  facet_wrap(tool~.) +
  labs(
    title = "Distribution of True Positives and False Positives per Gene and Tool",
    x = "Tool",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()

joined_tp <- left_join(joined_tp, gene_pp, by=c("gene_id"))
joined_tp <- joined_tp %>% filter(PP!=0)

joined_tp <- joined_tp %>%
  mutate(category = case_when(
    TP == PP & TP!=0 & FP == 0 ~ "All detected",
    TP > 0 & FP == 0 ~ "Only matching TP",
    TP > 0 & FP > 0  ~ "TP and FP",
    TP == 0 & FP > 0 ~ "Only FP",
    TRUE             ~ "None"   # Covers the case where TP == 0 and FP == 0
  ))

# View the updated data frame

png(sprintf("./tx_noise_fig/tx_gene_cat.png"), width=5000, height=2000, res=300)
ggplot(joined_tp %>% filter(category!="None"), aes(x = tool, fill = category)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Distribution of Gene Categories by Tool",
    x = "Tool",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



tx0 <- truthfiles_tx <- read.csv(paste0(sprintf(outdir, "pair", "4", "r1", "_0"), "/results/salmon_res_tx_N_T.txt"), sep="\t")
tx1 <- truthfiles_tx <- read.csv(paste0(sprintf(outdir, "pair", "8", "r1", "_0.1"), "/results/truthtable_tx.csv"), sep="\t")
tx0 %>% group_by(gene_id) %>% summarise(status=sum(status), events=max(events)) %>% group_by(events) %>% summarise(sum(status)) 
