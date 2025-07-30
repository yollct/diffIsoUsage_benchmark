library(tidyverse)


res_data <- read_delim("/nfs/scratch/chit/new_simulations/pair_50_4_r1_0/results/rsem_res_tx_N_T.txt")  # P-Werte auslesen

truth_data <- read_delim("/nfs/scratch/chit/new_simulations/pair_50_4_r1_0/results/truthtable_tx.csv", delim = "\t")

dim(res_data)
dim(truth_data)

sim = "kal"

# Function to calculate Precision and Recall
calculate_metrics <- function(res_data, truth_data, tool_column) {
 
  # Extract p values for the specific tools
  res_data_filtered <- res_data %>%
    filter(!!sym(tool_column) < 0.05)  # Signifikante Treffer für das Tool
  print(dim(res_data_filtered))
  # True Positives: Number of genes, which are labeled as positive in truth_data and occur also in res_data_filtered
  true_positives <- sum(truth_data$status[truth_data$feature_id %in% res_data_filtered$feature_id] == 1)
  
  # False Positives: Number of genes labeled as negative in truth_data and present in res_data_filtered
  false_positives <- sum(truth_data$status[truth_data$feature_id %in% res_data_filtered$feature_id] == 0)
  
  # False Negatives: Number of genes labeled as positive in truth_data but not present in res_data_filtered
  false_negatives <- sum(truth_data$status == 1 & !truth_data$feature_id %in% res_data_filtered$feature_id)
  
  precision <- ifelse((true_positives + false_positives) > 0, 
                      true_positives / (true_positives + false_positives), 0)
  recall <- ifelse((true_positives + false_negatives) > 0, 
                   true_positives / (true_positives + false_negatives), 0)
  
  return(c(precision = precision, recall = recall))
}


process_condition <- function(condition_folder) {
  tools <- c("iso_ktsp", "drimseq", "dexseq", "dturtle", "seqGSEA", "cuffdiff", "junctionseq",
             "saturn", "DSGseq", "nbsplice", "LimmaDS", "edgeR")
  

  res_file <- file.path(condition_folder, "results", sprintf("%s_res_tx_N_T.txt", sim))
  truth_file <- file.path(condition_folder, "results", "truthtable_tx.csv")
  

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

df
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
  geom_point(position = position_dodge(width = 0.4), size = 1) +
  geom_errorbar(aes(ymin = mean_precision - sd_precision, ymax = mean_precision + sd_precision),
                position = position_dodge(width = 0.4), width = 0.2) +
  labs(title = "Mean Precision by Noise Level and Condition",
       x = "Tool",
       y = "Mean Precision") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold"),
        strip.text.y = element_text(size = 10, face = "bold")) +
  facet_grid(rows = vars(condition_base), cols = vars(noise_level), scales = "free_y")

# Print the plot
print(precision_plot)

# Save the recall plot to a file
ggsave(sprintf("/nfs/proj/is_benchmark/tx_noise_fig/meanper_%s.png", sim), 
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
  geom_point(position = position_dodge(width = 0.4), size = 1) +
  geom_errorbar(aes(ymin = mean_recall - sd_recall, ymax = mean_recall + sd_recall),
                position = position_dodge(width = 0.4), width = 0.2) +
  labs(title = "Mean Recall by Noise Level and Condition",
       x = "Tool",
       y = "Mean Recall") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 10, face = "bold"),
        strip.text.y = element_text(size = 10, face = "bold")) +
  facet_grid(rows = vars(condition_base), cols = vars(noise_level), scales = "free_y")

# Print the plot
print(recall_plot)

# Save the recall plot to a file
ggsave(sprintf("/nfs/proj/is_benchmark/tx_noise_fig/meanre_%s.png", sim), 
       plot = last_plot(),         # Saves the last plot generated
       width = 8,                  # Set the width of the plot
       height = 6,                 # Set the height of the plot
       dpi = 300,
       bg = "white")   
