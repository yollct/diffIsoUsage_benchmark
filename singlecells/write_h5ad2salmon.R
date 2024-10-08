library(tidyverse)

# Assuming count_matrix is your single-cell count matrix
# Replace this with your actual matrix
count_matrix <- matrix(rnorm(100), nrow=10)  # Example data
rownames(count_matrix) <- paste0("Gene", 1:nrow(count_matrix))  # Example gene names
colnames(count_matrix) <- paste0("Sample", 1:ncol(count_matrix))  # Example sample names

new_count <- readRDS("/nfs/home/students/chit/is_benchmark/singlecells/rawdata/simulation_scd/sim_data_balance/sim_balance_2ct_rep50.rds")

# Function to write each column of the count matrix to a Salmon-like quant.sf file
write_salmon_quant_files <- function(count_matrix) {
  for (sample in colnames(count_matrix)) {
    # Extract counts for the current sample
    counts <- count_matrix[, sample]
    
    # Create a data frame similar to Salmon's quant.sf
    quant_sf <- tibble(
      Name = rownames(count_matrix),
      Length = NA,  # Assuming Length is not available
      EffectiveLength = NA,  # Assuming EffectiveLength is not available
      TPM = NA,  # Assuming TPM is not calculated here
      NumReads = counts
    )
    dir.create("test")
    dir.create(paste0("test/",sample))
    # File path for the output
    file_path <- paste0("test/", sample, "/quant.sf")
    
    # Write the data frame to a file
    write_tsv(quant_sf, file_path, col_names = TRUE)
  }
}

# Run the function with your count matrix
write_salmon_quant_files(count_matrix)
