# compute the similarity matrix using the euclidean distance
# and the affinity matrix using the Gaussian kernel.

# Load libraries
library(proxy)
library(SNFtool)
library(data.table)
library(pheatmap)

# 1. Data Preparation
file_path <- "subset_tables/subset_table.tsv"

# Read the TSV file, skipping the first descriptive title row
# - header = TRUE: the first row after skipping is treated as header
cat("Reading data...\n")
mydata <- fread(file_path, header = TRUE, sep = "\t", skip = 1)
print(dim(mydata))  # Expected: 39 rows x 9512 columns

# Extract OTU IDs and OTU counts
otu_ids <- mydata[[1]]  # First column: OTU IDs
otu_counts <- as.matrix(mydata[, -1])  # Remove the first column to get OTU counts # nolint
rownames(otu_counts) <- otu_ids  # Assign OTU IDs as row names
sample_data <- t(otu_counts)  # Dimensions: 9511 samples x 39 OTUs
sample_dataframe <- as.data.frame(sample_data)
sample_datarame <- apply(sample_dataframe, 2, as.numeric)
print(dim(sample_dataframe))  # Expected: 9511 x 39

# Visualize the uclidean distance matrix
subset_size <- 100
sample_dataframe_subset <- sample_dataframe[1:subset_size, ]
dist_matrix <- as.matrix(dist(sample_dataframe_subset, method = "euclidean"))

# transform Distances to Similarity Scores
sigma <- 0.5
similarity_matrix <- exp(-dist_matrix^2 / (2 * sigma^2))



# # SMF matrix
# W <- affinityMatrix(dist_matrix, K = 20, sigma = 0.5)
# dim(W)  # Expected: subset_size x subset_size

























# Cosine
# -----------------------------------------------
# Final R Script for Computing Patient Similarity Matrix
# -----------------------------------------------

# Load libraries
library(proxy)
library(SNFtool)
library(data.table)
library(Matrix)
library(ggplot2)
library(pheatmap)
library(umap)
library(doParallel)
library(foreach)

# 1. Data Preparation
file_path <- "subset_tables/subset_table.tsv"

# Read the TSV file, skipping the first descriptive title row
# - header = TRUE: the first row after skipping is treated as header (sample names)
# - sep = "\t": tab-separated values
# - skip = 1: skip the first descriptive row
cat("Reading data...\n")
mydata <- fread(file_path, header = TRUE, sep = "\t", skip = 1)
print(dim(mydata))  # Expected: 39 rows x 9512 columns

# Extract OTU IDs and OTU counts
otu_ids <- mydata[[1]]  # First column: OTU IDs
otu_counts <- as.matrix(mydata[, -1])  # Remove the first column to get OTU counts
rownames(otu_counts) <- otu_ids  # Assign OTU IDs as row names

# Transpose the data to have samples as rows and OTUs as columns
sample_data <- t(otu_counts)  # Dimensions: 9511 samples x 39 OTUs

# Convert to numeric
cat("Converting data to numeric...\n")
sample_data <- apply(sample_data, 2, as.numeric)

# Verify the dimensions
print(dim(sample_data))  # Expected: 9511 x 39

# Check for NA values
if (any(is.na(sample_data))) {
  stop("NA values found in sample data. Please check the input file for non-numeric entries.")
}

# 2. Creating Similarity Networks

# Since there are only 39 OTUs, compute a single similarity matrix without splitting into groups
# Normalize data if necessary (e.g., z-score normalization)
cat("Normalizing data...\n")
sample_data_norm <- scale(sample_data)

# Handle any potential NA values after scaling
if (any(is.na(sample_data_norm))) {
  stop("NA values found after normalization. Check for constant columns.")
}

# 3. Compute Cosine Similarity Matrix

cat("Computing cosine similarity matrix...\n")
similarity_matrix <- simil(sample_data_norm, method = "cosine", diag = TRUE, upper = TRUE)
similarity_matrix <- as.matrix(similarity_matrix)
print(dim(similarity_matrix))  # Should be 9511 x 9511

cat("Generating heatmap for the similarity matrix...\n")
pheatmap(similarity_matrix,
         cluster_rows = TRUE,   # Cluster rows to group similar samples
         cluster_cols = TRUE,   # Cluster columns similarly
         show_rownames = TRUE,  # Show sample names on rows
         show_colnames = TRUE,  # Show sample names on columns
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste("Cosine Similarity Matrix (", subset_size, "Samples)", sep = ""))

# Save the similarity matrix
saveRDS(similarity_matrix, "cosine_similarity_matrix.rds")
cat("Cosine similarity matrix saved as 'cosine_similarity_matrix.rds'\n")

# 4. Clustering and Reordering Samples

# a. Spectral Clustering using SNFtool
k <- 3  # Number of clusters
cat("Performing spectral clustering with", k, "clusters...\n")
clusters <- spectralClustering(similarity_matrix, K = k)

# Create a data frame with cluster assignments
clustered_samples <- data.frame(Sample = rownames(similarity_matrix), Cluster = clusters)

# View the first few cluster assignments
head(clustered_samples)

# Save clustering results to disk
write.csv(clustered_samples, "clustered_samples.csv", row.names = FALSE)
cat("Cluster assignments saved as 'clustered_samples.csv'\n")

# b. Reordering the Similarity Matrix based on clusters
cat("Reordering the similarity matrix based on cluster assignments...\n")
ordered_indices <- order(clusters)
ordered_similarity <- similarity_matrix[ordered_indices, ordered_indices]
print(dim(ordered_similarity))  # Should still be 9511 x 9511

# Save the ordered similarity matrix (optional)
saveRDS(ordered_similarity, "ordered_similarity_matrix.rds")
cat("Ordered similarity matrix saved as 'ordered_similarity_matrix.rds'\n")

# 5. Visualizing the Similarity Matrix

# a. Sample Heatmap
subset_size <- 1000  # Adjust based on system capabilities
if (subset_size > nrow(ordered_similarity)) {
  stop("Subset size exceeds the number of available samples.")
}

subset_indices <- ordered_indices[1:subset_size]
subset_similarity <- ordered_similarity[subset_indices, subset_indices]

cat("Generating heatmap for a subset of samples...\n")
pheatmap(subset_similarity,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = paste("Cosine Similarity Matrix (Subset of", subset_size, "Samples)"))

# b. Dimensionality Reduction with UMAP

cat("Applying UMAP on the subset distance matrix...\n")
distance_matrix_subset <- 1 - subset_similarity  # Convert similarity to distance

umap_config <- umap.defaults
umap_config$n_neighbors <- 15
umap_config$min_dist <- 0.1

umap_result <- umap(as.dist(distance_matrix_subset), config = umap_config)

# Plot UMAP results with cluster coloring
plot(umap_result$layout, 
     col = clusters[subset_indices],
     pch = 19, 
     main = paste("UMAP Projection of", subset_size, "Samples"),
     xlab = "UMAP 1", 
     ylab = "UMAP 2",
     col.main = "black")

# Add a legend (optional)
legend("topright", legend = paste("Cluster", 1:k), 
       col = 1:k, pch = 19, title = "Clusters")
