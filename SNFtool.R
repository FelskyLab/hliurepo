# Load required libraries
library(data.table)
library(readr)
library(tibble)
library(dplyr)
library(SNFtool)

taxonomy <- 6
file_name_vector <- c("./subset_tables/subset_table_phylum.tsv",
                      "./subset_tables/subset_table_class.tsv",
                      "./subset_tables/subset_table_order.tsv",
                      "./subset_tables/subset_table_family.tsv",
                      "./subset_tables/subset_table_genus.tsv",
                      "./subset_tables/subset_table_species.tsv")

# Initialize a list to store W matrices
W_list <- list()

for (i in 1:taxonomy) {
  file_path <- file_name_vector[i]

  # 1. Data Preparation
  mydata <- fread(file_path, header = TRUE, sep = "\t", skip = 1)
  print(dim(mydata))  # Expected: 39 rows x 9512 columns

  otu_ids <- mydata[[1]]
  otu_counts <- as.matrix(mydata[, -1])
  rownames(otu_counts) <- otu_ids
  sample_data <- t(otu_counts)
  sample_dataframe <- as.data.frame(sample_data)

  # 2. SNF Analysis
  subset_size <- 100
  sample_dataframe_subset <- sample_dataframe[1:subset_size, ]

  dist_matrix <- as.matrix(dist(sample_dataframe_subset))
  W <- affinityMatrix(dist_matrix, K = 20, sigma = 0.5)

  # Store the W matrix in the list
  W_list[[i]] <- W
}

# Assign W1 to W6 from the list
W1 <- W_list[[1]]
W2 <- W_list[[2]]
W3 <- W_list[[3]]
W4 <- W_list[[4]]
W5 <- W_list[[5]]
W6 <- W_list[[6]]

# Perform SNF on the list of W matrices
W <- SNF(list(W1, W2, W3, W4, W5, W6), K = 20, t = 20)

# open a PDF device
pdf("clustering_heatmap.pdf", width = 8, height = 6)

# Perform spectral clustering
clusters <- spectralClustering(W, K = 3)

# Display clusters with heatmap
dummy <- displayClustersWithHeatmap(W, clusters)

# Close the PDF device
dummy <- dev.off()