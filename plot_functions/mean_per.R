library(tidyverse)

source("/nfs/proj/is_benchmark/plot_functions/new_functions.R")

res_data <- read_delim("/nfs/scratch/chit/new_simulations/pair_50_4_r1_0/results/rsem_res_gene_N_T.txt")  # P-Werte auslesen

truth_data <- read_delim("/nfs/scratch/chit/new_simulations/pair_50_4_r1_0/results/truthtable_gene.csv", delim = "\t")

dim(res_data)
dim(truth_data)

sim = "kal"

thresholds <- data.frame(tools=c("iso_ktsp","drimseq","dexseq","dturtle","seqGSEA","cuffdiff","junctionseq","saturn","drimseq_stageR", "dexseq_stageR","saturn_stageR", "DSGseq", "nbsplice", "LimmaDS", "edgeR", "BANDITS"),
                        thres=c(0.5,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05, 5, 0.05, 0.05, 0.05, 0.05))
row.names(thresholds) <- thresholds$tools

# Function to calculate Precision and Recall
calculate_metrics <- function(methoddf, truthfile, method) {
 
  methoddf <- methoddf[methoddf$feature_id %in% truthfile$feature_id,]
  namet <- methoddf[,method]
  
  names(namet) <- "value"
  namet$feature_id <- methoddf$feature_id
  x<-thresholds[method,]$thres
  

  namet$value <- lapply(namet$value, function(y){ifelse(x==0.05, ifelse(is.na(y), 1, y), ifelse(is.na(y), 0, y))}) %>% unlist()
  
  #precisions <- sapply(thresholds, function(x){
  boomet <- lapply(namet$value, function(y){ifelse(x==0.05, y<x, y>x)}) %>% unlist

  tp <- sum(unique(namet$feature_id[boomet]) %in% truthfile[truthfile$status==1,]$feature_id) 
  fp <- sum(!unique(namet$feature_id[boomet]) %in% truthfile[truthfile$status==1,]$feature_id) 
  fn <- sum(!unique(namet$feature_id[!boomet]) %in% truthfile[truthfile$status==0,]$feature_id) 

  pp <- length(unique(namet$feature_id[boomet]))
  precision <- tp/pp

    #   })
  #print(truthfile$feature_id[!truthfile[truthfile$status==1,]$feature_id %in% names(met)[met<x]])
  #recall <- sapply(thresholds, function(x){
  
  p <- length(truthfile[truthfile$status==1,]$feature_id)
  recall <- tp/p

  f1 <- 2*tp/(2*tp+fp+fn)
  
  fdr <- fp / (tp + fp) # False Discovery Rate calculation
  
  # outout <- data.frame(
  #   metric_value = c(precision, recall, f1, fdr, tp, pp, length(truthfile[truthfile$status == 1,]$feature_id)),
  #   metric_type = c("precision", "recall", "f1", "fdr", "true_pos", "detected_pos", "real_pos"),
  #   tool = rep(method, 7),
  #   thresholds = rep(x, 7)
  # )
  #    return(tp/p)})


  return(c(precision = precision, recall = recall))
}


process_condition <- function(condition_folder) {
  tools <- c("iso_ktsp", "drimseq", "dexseq", "dturtle", "seqGSEA", "cuffdiff", "junctionseq",
             "saturn", "DSGseq", "nbsplice", "LimmaDS", "edgeR", "BANDITS")
  

  res_file <- file.path(condition_folder, "results", sprintf("%s_res_gene_N_T.txt", sim))
  truth_file <- file.path(condition_folder, "results", "truthtable_gene.csv")
  

  if (!file.exists(res_file) || !file.exists(truth_file)) {
    return(NULL)
  }
  

  res_data <- read_delim(res_file)
  truth_data <- read_delim(truth_file, delim = "\t")
  

  metrics <- map_dfr(tools, function(tool) {
    tool_column <- tool  # name of tool is where the p value lies
    print(tool_column)
    print("next_round")
    if (tool_column %in% colnames(res_data)) {
      metrics <- calculate_metrics(res_data, truth_data, tool_column)
      return(data.frame(tool = tool, precision = metrics["precision"], recall = metrics["recall"]))
    } else {
      return(data.frame(tool = tool, precision = NA, recall = NA)) 
    }
  })
  
  metrics$condition <- basename(condition_folder)
  return(metrics)
}

process_all_conditions <- function(main_folder) {
  condition_folders <- list.dirs(main_folder, recursive = FALSE)  
  
  all_metrics <- map_dfr(condition_folders, process_condition)
  return(all_metrics)
}

main_folder <- "/nfs/scratch/chit/new_simulations"  
all_data <- process_all_conditions(main_folder)

df <- all_data

df <- df %>% mutate(condition_base = sub("_(r[1-5])_", "_", condition))

df$condition_base <- sub("_\\d+(\\.\\d+)?$", "", df$condition_base)

df$noise_level <- sub(".*_(.*)$", "\\1", df$condition)

df <- df %>% filter(!is.na(precision) & !is.na(recall))


# Load required libraries
library(ggplot2)
library(dplyr)

# Prepare the data by calculating mean precision for each tool and condition combination
precision_summary <- df %>%
  group_by(condition_base, tool, noise_level) %>%
  summarise(mean_precision = mean(precision),
            sd_precision = sd(precision), .groups = "drop")

# Create the plot with faceting on noise levels (5 columns) and major conditions (4 rows)
precision_plot <- ggplot(precision_summary, aes(x = tool, y = mean_precision, color = tool)) +
  geom_point(position = position_dodge(width = 0.4), size = 0.6) +
  geom_errorbar(aes(ymin = mean_precision - sd_precision, ymax = mean_precision + sd_precision),
               position = position_dodge(width = 0.4), width = 0.2) +
  labs(title = "Mean Precision by Noise Level and Condition",
       x = "Tool",
       y = "Mean Precision") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        strip.text.x = element_text(size = 10, face = "bold"),
        strip.text.y = element_text(size = 10, face = "bold")) +
  facet_grid(rows = vars(condition_base), cols = vars(noise_level), scales = "free_y")

# Print the plot
print(precision_plot)

# Save the recall plot to a file
ggsave(sprintf("/nfs/proj/is_benchmark/noise_fig/meanper_%s.png", sim), 
       plot = last_plot(),         # Saves the last plot generated
       width = 8,                  # Set the width of the plot
       height = 6,                 # Set the height of the plot
       dpi = 300,
       bg = "white")   

# Prepare the data by calculating mean recall for each tool and condition combination
recall_summary <- df %>%
  group_by(condition_base, tool, noise_level) %>%
  summarise(mean_recall = mean(recall),
            sd_recall = sd(recall), .groups = "drop")

# Create the plot with faceting on noise levels (5 columns) and major conditions (4 rows)
recall_plot <- ggplot(recall_summary, aes(x = tool, y = mean_recall, color = tool)) +
  geom_point(position = position_dodge(width = 0.4), size = 0.6) +
  geom_errorbar(aes(ymin = mean_recall - sd_recall, ymax = mean_recall + sd_recall),
                position = position_dodge(width = 0.4), width = 0.2) +
  labs(title = "Mean Recall by Noise Level and Condition",
       x = "Tool",
       y = "Mean Recall") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        strip.text.x = element_text(size = 10, face = "bold"),
        strip.text.y = element_text(size = 10, face = "bold")) +
  facet_grid(rows = vars(condition_base), cols = vars(noise_level), scales = "free_y")

# Print the plot
print(recall_plot)

# Save the recall plot to a file
ggsave(sprintf("/nfs/proj/is_benchmark/noise_fig/meanre_%s.png", sim), 
       plot = last_plot(),         # Saves the last plot generated
       width = 8,                  # Set the width of the plot
       height = 6,                 # Set the height of the plot
       dpi = 300,
       bg = "white")   
