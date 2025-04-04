library(data.table)
library(ggplot2)
library(compositions)

# the load_data function will only accept one argument as input. It then
# does a similar dataprocessing step so it will return a dataframe. The
# rows are the samples and the columns are the OTUs.
load_data_stacked <- function(file_path) {
  # Step 1: Data Preparation
  # Note fread is faster than readr::read_tsv
  mydata <- fread(file_path, header = TRUE, sep = "\t", skip = 1)
  # Examine the original data. There are 9511 samples corresponding
  # to 9512 columns (1 column for the OTU ID and 9511 columns for the samples)
  cat("Original data dimensions: ", dim(mydata), "\n")
  # Expected: 39 rows x 9512 columns

  # The first colomn is the OTU ID (bacteria kingdom/class/...)
  otu_ids <- mydata[[1]]
  otu_counts <- as.matrix(mydata[, -1])
  rownames(otu_counts) <- otu_ids
  # We do a matrix transpose so the sample IDs are the rows now
  sample_data <- t(otu_counts)
  # We store the result as a dataframe
  sample_dataframe <- as.data.frame(sample_data)
  return(sample_dataframe)
}

sample_dataframe <- load_data_stacked("subset_tables/subset_table_order.tsv")
subset_size <- 1000
subset_df <- sample_dataframe[1:subset_size, ]

# Optional: Apply CLR transformation for compositional data
subset_clr <- clr(as.matrix(subset_df))

# Perform PCA
pca_result <- prcomp(subset_df, center = TRUE, scale. = FALSE)
# For CLR-transformed data:
# pca_result <- prcomp(subset_clr, center = TRUE, scale. = FALSE)

# Summary of PCA
summary(pca_result)

# Extract PCA scores
pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$SampleID <- rownames(subset_df)

# Basic PCA plot
p <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(color = "blue", alpha = 0.6) +
  labs(
    title = "PCA of Bacterial Relative Abundance",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "% Variance)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "% Variance)")
  ) +
  theme_minimal()

ggsave("PCA.pdf",
       plot = p,
       width = 15,
       height = 10)
