library(tidyverse)
library(ggplot2)
library(UpSetR)
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(GGally)
library(fmsb)
library(gridExtra)
source("/home/chit/norm_is/plot_functions/new_functions.R")

path <- "/home/chit/norm_is/tx_noise_fig"
outdir <- "/home/chit/norm_chit/new_simulations/%s_50_%s_%s%s"

noises <- c("_0", "_0.1","_0.5", "_0.7", "_0.9")

thisrep <- c('4','8')
seqtype <- 'pair'
whichrep <- 'r1'

thresholds <- data.frame(tools=c("iso_ktsp","drimseq","dexseq","dturtle","seqGSEA","cuffdiff","junctionseq","saturn","drimseq_stageR", "dexseq_stageR","saturn_stageR", "DSGseq", "nbsplice", "LimmaDS", "edgeR", "BANDITS"),
                        thres=c(0.5,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05, 5, 0.05, 0.05, 0.05, 0.05))
row.names(thresholds) <- thresholds$tools

get_overall <- function(split){
    overalldf <- do.call(rbind, lapply(whichrep, function(R){
            do.call(rbind, lapply(seqtype, function(seq){
            do.call(rbind, lapply(thisrep, function(rep){
                do.call(rbind, lapply(noises, function(noise){
                    print(paste0(sprintf(outdir, seq, rep, R, noise)))
                    #print(paste0(outdir, seq, rep, R, noise))
                    #print(seq)
                    truthfiles_g <- read.csv(paste0(sprintf(outdir, seq, rep, R, noise), "/results/truthtable_gene.csv"), sep="\t") ### change to result
                    truthfiles_tx <- read.csv(paste0(sprintf(outdir, seq, rep, R, noise), "/results/truthtable_tx.csv"), sep="\t")
                    
                    
                    if (is.null(split)){
                        print("not consider DTE")
                        # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                    } else {
                        if (split!="events"){
                        print("not consider DTE")
                        # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status=ifelse(events=="DTE", 0, status))
                        }
                    }

                    diffiso <- read.csv(paste0(sprintf(outdir, seq, rep, R, noise), "/results/isoinfo.csv"), sep="\t")
                    isorat <- diffiso %>% group_by(gene_id) %>% summarise(diffiso=max(diffiso))
                    truthfiles_g$diffiso <- isorat$diffiso
                    truthfiles_g <- truthfiles_g %>% dplyr::mutate(diffiso=ifelse(diffiso>0.8, ">0.8", ifelse(diffiso>=0.6, "0.6-0.8", ifelse(diffiso>=0.5, "0.5-0.6", ifelse(diffiso>=0.4, "0.4-0.5", ifelse(diffiso>=0.3, "0.3-0.4", ifelse(diffiso>=0.2, "0.2-0.3", ifelse(diffiso>=0.1, "0.1-0.2", "0-0.1")))))))) 
                    print(truthfiles_g %>% nrow)
                    # truthfiles_g <- truthfiles_g %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    # truthfiles_tx <- truthfiles_tx %>% dplyr::mutate(status = ifelse(events=="DTE", 0, status))
                    
                    do.call(rbind, lapply(c("salmon", "kal", "rsem"), function(quant){
                        print(paste0(sprintf(outdir, seq,rep, R, noise), "/results/", sprintf("%s_res_tx_N_T.txt", quant)))
                        df_g <- read_delim(paste0(sprintf(outdir, seq,rep, R, noise), "/results/", sprintf("%s_res_tx_N_T.txt", quant)))
                        print( colnames(df_g))
                        thisthres <- thresholds[colnames(df_g)[2:ncol(df_g)],]$thres
                        df_g <- df_g %>% group_by(feature_id) %>% summarise(across(colnames(select(df_g, -feature_id)), .fns = min), .groups = "keep")

                        df_g <- df_g[grepl("ENST", df_g$feature_id),]

                        outputpr <- cal_pre_re(df_g, truthfiles_tx, split=split)
                        
                        this_pr <- pivot_output(outputpr, split=split)
                        this_pr$quant_tool <- quant
                        this_pr$rep <- paste(as.character(rep), 'replicates')
                        this_pr$seq <- seq
                        this_pr$sim <- R
                        this_pr$con <- paste0(quant, "_", paste(as.character(rep), 'replicates'))
                        print(noise)
                        this_pr$noise <- ifelse(noise == "_0", 0, 
                                                ifelse(noise == "_0.1", 0.1, 
                                                       ifelse(noise == "_0.5", 0.5, 
                                                              ifelse(noise == "_0.7", 0.7, 
                                                                     ifelse(noise == "_0.9", 0.9, NA))))) 
                        this_pr
                    }))
                }))
            }))
        }))
    }))
    overalldf$fdr <- lapply(overalldf$fdr, function(x){as.numeric(x)}) %>% unlist
    overalldf$f1 <- lapply(overalldf$f1, function(x){as.numeric(x)}) %>% unlist
    overalldf
}

#diffiso <- read.csv("/home/rstudio/benchmark_paper/sample_result/pair_50_4_r1_0/results/isoinfo.csv", sep="\t")
#diffiso <- read_delim("/home/rstudio/benchmark_paper/sample_result/pair_50_4_r1_0/results/salmon_res_gene_N_T.txt")


colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9), "#008000" , "black" )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4),"#008000" , "black" )

# where will this be used?
seq_label <- as_labeller(c('pair'="Paired-end", 'single'='Single-end', '2'='2', '3'='3','4'='4','5'='5', '2-4'='2-4','5-9'='5-9','>9'='>9','DTE'='S1','DTU'='S2','IS'='S3','0-0.1'='0-0.1', '0.1-0.2'='0.1-0.2','0.3-0.4'='0.3-0.4', '0.4-0.5'='0.4-0.5', '0.5-0.6'='0.5-0.6', '0.6-0.8'='0.6-0.8', '>0.8'='>0.8', '8 replicates'='8 replicates', '4 replicates'='4 replicates', '0'='0','0.5'='0.5', '0.1'='0.1', "kal" = "Kallisto", "rsem" = "RSEM","salmon" = "Salmon", 'kal_8 replicates'='Kallisto: 8 replicates', 'kal_4 replicates'='Kallisto: 4 replicates', 'rsem_8 replicates'='RSEM: 8 replicates', 'rsem_4 replicates'='RSEM: 4 replicates', 'salmon_8 replicates'='Salmon: 8 replicates', 'salmon_4 replicates'='Salmon: 4 replicates'))

alldf <- get_overall(split=NULL)


overalldffdr <- alldf %>% dplyr::filter(sim==whichrep&seq==seqtype) 
overalldffdr_no_na <- overalldffdr[!is.na(overalldffdr$f1), ]

#write.table(f1score, "/nfs/scratch/chit/GSE222260/analysis/f1scores.csv", sep="\t")

plot_list <- c()
plot_name <- c()
n<-1

for (rps in c("4 replicates", "8 replicates")) {
  for (spl in sort(unique(overalldffdr_no_na$quant_tool))) {
    data <- overalldffdr_no_na %>%
      filter(quant_tool == spl & rep == rps) %>%
      select(noise, f1, tool) %>%
      pivot_wider(names_from = tool, values_from = f1, id_cols = noise, values_fn = mean) %>%
      as.data.frame()
    
    row.names(data) <- paste0("bg-", data$noise)
    print(data)
    data <- rbind(rep(0.6, ncol(data)), rep(0, ncol(data)), data)
    data <- data %>% select(-noise)
    
    # Ensure all columns are numeric
    data[] <- lapply(data, as.numeric)
    
    plot_list[[n]] <- data
    plot_name[[n]] <- paste0(rps, "_", spl)
    n <- n + 1
  }
}


newdata <- overalldffdr_no_na %>%
      select(quant_tool, rep, con, noise, f1, tool) %>%
      as.data.frame() %>% 
      mutate(tool=fct_reorder(tool, f1))
# Compute mean value per category


# Reorder factor levels by mean
# For each Object, sort Category by mean value




png(sprintf("%s/alltools_allf1radar_%s_%s.png", path, seqtype, whichrep), width=4000, height=2000, res=300)
newdata$noise <- as.character(newdata$noise)
# ggplot(newdata, aes(x=tool, y=f1, shape=noise))+
#   geom_point(size=5, alpha=0.8)+geom_path(aes(x=tool, y=f1, group=tool))+
#   facet_grid(quant_tool ~ rep, labeller=seq_label)+
#   theme_minimal(base_size = 14) +
#   scale_color_manual(values=colMap)+
#   theme(axis.text=element_text(size=20), axis.text.x =  element_blank(), axis.title = element_text(size=20),
#         strip.text=element_text(size=20), legend.text=element_text(size=15), legend.title=element_text(size=15)) +
#         labs(color="Tools", shape="Background")+
#         ggtitle("F1 in different event scenarios in paired-end data")
# dev.off()

nrows<-2
ncols<-3
par(mfrow = c(nrows, ncols), mar = c(0.5,0.5,0.5,0.5), oma=c(4,4,4,4)) # 4 rows, 2 columns

for (i in seq_along(plot_list)) {
radarchart(plot_list[[i]], axistype=1 ,
        #custom polygon
        pcol=colors_border ,  plwd=4 , plty=1,
        #custom the grid
        cglcol="grey", cglty=1, axislabcol="#5E5E5E", caxislabels=seq(0,0.6,0.15), cglwd=0.8,
        #custom labels
        vlcex=1.5, calcex=1

 )
}
#Colors/shades for each group (5 rows)
# pch_list <- c(16, 17, 15, 18, 8)  # point shapes
# col_list <- gray.colors(5, start = 0.2, end = 0.8)  # colors for lines/points

# for (i in seq_along(plot_list)) {
#   data_i <- plot_list[[i]]  # 5 rows (groups) × N columns (categories)
#   num_groups <- nrow(data_i)
#   num_cats <- ncol(data_i)
#   category_labels <- colnames(data_i)

#   # Set up empty plot space
#   plot(1, type = "n",
#        xlim = c(1, num_cats),
#        ylim = c(0, 1),
#        xaxt = "n",
#        yaxt = TRUE,
#        xlab = "",
#        ylab = "",
#        main = "")

#   axis(1, at = 1:num_cats, labels = category_labels, las = 2, cex.axis = 1.2)
#   axis(2, las = 1, cex.axis = 1.2)
#   box()

#   # Plot each group as a path
#   for (g in 1:num_groups) {
#     values <- as.numeric(data_i[g, ])
#     points(1:num_cats, values, type = "b", lwd = 2,
#            pch = pch_list[g], col = col_list[g])
#   }
# }

# legend("bottom", legend = paste("Group", 1:5),
#        pch = pch_list, pt.cex = 1.5, horiz = TRUE, bty = "n", inset = -0.2, xpd = NA, cex = 1.2)

# str(plot_list[[i]])



row_labels <- c("4 replicates", "8 replicates") # Adjust as needed
for (i in 1:nrows) {
  mtext(row_labels[i], side = 2, at = (nrows - i + 0.5)/nrows, srt=90, outer = TRUE, line = 1, cex=1.8)
}

# Add column labels
col_labels <- c("Kallisto", "RSEM", "Salmon")  # Adjust as needed
for (i in 1:ncols) {
  mtext(col_labels[i], side = 3, at = (i - 0.5)/ncols, las = 1, outer = TRUE, line = 1, cex=1.8)
}

# 
dev.off()

library(ggplot2)
library(dplyr)

# Calculate mean FDR and mean recall for each quant_tool and rep, averaged over simulation replicates
fdr_recall_summary <- overalldffdr_no_na %>%
  group_by(quant_tool, rep, tool, noise) %>%
  summarise(mean_recall = mean(recall), 
            mean_fdr = mean(fdr), .groups = "drop")

library(ggplot2)
library(dplyr)

# Assuming your tibble is named fdr_recall_summary
# Create the plot
ggplot(fdr_recall_summary, aes(x = mean_recall, y = mean_fdr, color = tool, shape = factor(noise))) +
  geom_point(size = 3) +   # Add points for each noise level
  geom_line(aes(group = interaction(tool)), size = 1) +    # Connect points by noise level and tool
  facet_grid(rep ~ quant_tool, scales = "free", space = "free") + 
  theme_minimal() +
  labs(
    x = "Recall",
    y = "FDR",
    title = "FDR vs. Recall for Different Tools and Replicates"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    strip.text = element_text(size = 12),
    strip.background = element_blank()
  ) 



