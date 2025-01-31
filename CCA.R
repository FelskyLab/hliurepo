# load libraries

library(data.table)
library(CCA)
library(CCP)
library(dplyr)

# 2) Define the custom function to read and transpose OTU table
# (same in stackedPlot.R)
load_data_stacked <- function(file_path) {
  # Use fread from data.table (fast reading).
  # The skip=1 below is because OTU table have 1 header line before columns.
  mydata <- fread(file_path, header = TRUE, sep = "\t", skip = 1)

  cat("Original data dimensions: ", dim(mydata), "\n")
  # e.g., "39 rows x 9512 columns" if you had 39 taxa
  # and 9511 samples + 1 for OTU ID
  # First column = OTU IDs
  otu_ids <- mydata[[1]]

  # The rest of the columns = sample counts/abundances
  otu_counts <- as.matrix(mydata[, -1])
  rownames(otu_counts) <- otu_ids
  # Transpose so rows = samples, columns = taxa
  sample_data <- t(otu_counts)

  # Convert to data frame
  sample_dataframe <- as.data.frame(sample_data)

  return(sample_dataframe)
}

# 3) filter the data
sample_dataframe <- load_data_stacked("subset_tables/subset_table_order.tsv")
threshold <- 0.001
col_means <- colMeans(sample_dataframe)
sample_dataframe_filtered <- sample_dataframe[, col_means > threshold, drop = FALSE]

# 4) load metadata and match sample IDs
metadata_file <- "dataset/metadata.tsv"
metadata <- fread(metadata_file, header = TRUE, sep = "\t") %>%
  as.data.frame()