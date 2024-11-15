# -----------------------------------------------
# Final R Script for Computing and Visualizing Similarity Matrix for First 50 Samples
# -----------------------------------------------


# Load libraries
library(proxy)
library(data.table)
library(Matrix)
library(ggplot2)
library(pheatmap)
library(umap)

# 1. Data Preparation
file_path <- "subset_tables/subset_table.tsv" # note this is the phylum level

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

# Select the first 50 samples
subset_size <- 200
# make srue subset_size is less than the number of samples
if (subset_size > nrow(sample_data)) {
  stop("Subset size exceeds the number of available samples.")
}

subset_data <- sample_data[1:subset_size, ]  # First 50 samples
cat("Selected the first", subset_size, "samples for analysis.\n")
print(dim(subset_data))  # Expected: 50 x 39

# Convert to numeric (if not already)
subset_data <- apply(subset_data, 2, as.numeric)

# Verify the dimensions
print(dim(subset_data))  # Should be 50 x 39

# 2. Computing the Similarity Matrix

# Normalize the data (optional but recommended for similarity measures)
# Z-score normalization
cat("Normalizing subset data...\n")
subset_data_norm <- scale(subset_data)

# Compute Cosine Similarity Matrix
cat("Computing cosine similarity matrix for the subset...\n")
similarity_matrix <- simil(subset_data_norm,
                           method = "cosine",
                           diag = TRUE,
                           upper = TRUE)

# Convert to a dense matrix
similarity_matrix <- as.matrix(similarity_matrix)

# Verify dimensions
print(dim(similarity_matrix))  # Should be 50 x 50

# 3. Visualizing the Similarity Matrix

# a. Heatmap Visualization
cat("Generating heatmap for the similarity matrix...\n")
pheatmap(similarity_matrix,
         cluster_rows = TRUE,   # Cluster rows to group similar samples
         cluster_cols = TRUE,   # Cluster columns similarly
         show_rownames = TRUE,  # Show sample names on rows
         show_colnames = TRUE,  # Show sample names on columns
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste("Cosine Similarity Matrix (", subset_size, "Samples)",
                      sep = ""))

# b. Dimensionality Reduction with UMAP (Optional)
# Convert similarity to distance
distance_matrix_subset <- 1 - similarity_matrix

cat("Applying UMAP on the subset distance matrix...\n")
umap_config <- umap.defaults
umap_config$n_neighbors <- 15
umap_config$min_dist <- 0.1

umap_result <- umap(as.dist(distance_matrix_subset), config = umap_config)






















# Plot UMAP results
cat("Plotting UMAP projection...\n")
plot(umap_result$layout, 
     col = "steelblue",
     pch = 19, 
     main = paste("UMAP Projection of", subset_size, "Samples"),
     xlab = "UMAP 1", 
     ylab = "UMAP 2",
     col.main = "black")

# Optionally, add sample labels
text(umap_result$layout, labels = rownames(subset_data), pos = 4, cex = 0.7)

# 4. Saving Results (Optional)

# Save the similarity matrix
saveRDS(similarity_matrix, "cosine_similarity_matrix_subset50.rds")
cat("Cosine similarity matrix for the subset saved as 'cosine_similarity_matrix_subset50.rds'\n")

# Save the heatmap (if not interactive)
# You can save the plot using the png or pdf device
png(filename = "cosine_similarity_heatmap_subset50.png", width = 800, height = 800)
pheatmap(similarity_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste("Cosine Similarity Matrix (", subset_size, "Samples)", sep = ""))
dev.off()
cat("Heatmap saved as 'cosine_similarity_heatmap_subset50.png'\n")

# Save the UMAP plot
png(filename = "UMAP_projection_subset50.png", width = 800, height = 600)
plot(umap_result$layout, 
     col = "steelblue",
     pch = 19, 
     main = paste("UMAP Projection of", subset_size, "Samples"),
     xlab = "UMAP 1", 
     ylab = "UMAP 2",
     col.main = "black")
text(umap_result$layout, labels = rownames(subset_data), pos = 4, cex = 0.7)
dev.off()
cat("UMAP projection saved as 'UMAP_projection_subset50.png'\n")
