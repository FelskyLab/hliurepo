#!/usr/bin/env Rscript
# 1) Load libraries
library(data.table)
library(CCA)
library(CCP)
library(dplyr)

# 2) Define your custom function to read and transpose OTU table
load_data_stacked <- function(file_path) {
  # Use fread from data.table (fast reading). Adjust skip if you have comment
  # lines.
  # The skip=1 below is because your OTU table might have 1 header line before
  # columns.
  mydata <- fread(file_path, header = TRUE, sep = "\t", skip = 1)
  cat("Original data dimensions: ", dim(mydata), "\n")
  # e.g., "39 rows x 9512 columns" if you had 39 taxa and
  # 9511 samples + 1 for OTU ID
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

# 3) Load the OTU table and filter low-abundance taxa
otu_file <- "subset_tables/subset_table_order.tsv"
# Load and transpose
sample_dataframe <- load_data_stacked(otu_file)
cat("After transpose, dimension is: ", dim(sample_dataframe), "\n")

# 3.5) Filter out taxa with average abundance < 1% (0.001 in fractional form)
# Filter out taxa with average abundance < 1% (0.01 in fractional form)
threshold <- 0.01  # 1%
col_means <- colMeans(sample_dataframe)
sample_dataframe_filt <- sample_dataframe[,
                                          col_means >= threshold,
                                          drop = FALSE]

cat("Dimension after filtering (< 1% avg): ", dim(sample_dataframe_filt), "\n")

# 4) Load metadata and exclude samples with missing BMI/age
metadata_file <- "dataset/metadata.tsv"
metadata <- fread(metadata_file, header = TRUE, sep = "\t") %>%
  as.data.frame()
# Set rownames in metadata to sample_name
rownames(metadata) <- metadata$sample_name
# Coerce bmi and age_years to numeric
# Any invalid string becomes NA automatically
metadata$bmi       <- as.numeric(as.character(metadata$bmi))
metadata$age_years <- as.numeric(as.character(metadata$age_years))
# EXCLUDE rows with NA in BMI or age_years 
# (this removes both originally NA and any invalid entries)
metadata_clean <- metadata %>%
  filter(!is.na(bmi) & !is.na(age_years))
cat("Metadata dimension after removing NA in BMI/age_years:", dim(metadata_clean), "\n")
# Now we want to match sample IDs in 'sample_dataframe_filt' to 'metadata_clean'
common_samples <- intersect(rownames(sample_dataframe_filt),
                            rownames(metadata_clean))
cat("Number of matching samples (after NA exclusion): ", length(common_samples), "\n")
# Subset both data frames to matching samples
otu_final <- sample_dataframe_filt[common_samples, , drop=FALSE]
metadata_final <- metadata_clean[common_samples, , drop=FALSE]


library(compositions)
library(ggplot2)

# Define an age group variable for grouping in PCA
metadata_final <- metadata_final %>%
  mutate(age_group = ifelse(age_years > 45, "Above45", "Below45"))

# Step A1: Basic zero replacement (pseudocount) - you could also use zCompositions
# If you have many zeros, consider a more refined approach (zCompositions, etc.)
pseudocount_value <- 1e-5
otu_versionA <- otu_final + pseudocount_value

# Step A2: Convert to 'acomp' class (compositions package expects that for clr())
otu_acomp <- acomp(otu_versionA)

# Step A3: CLR transform
clr_mat_A <- clr(otu_acomp)
# 'clr_mat_A' is now a matrix or data.frame of the CLR-transformed data

# Step A4: PCA
# Some users do prcomp(..., center=TRUE, scale.=FALSE) even after CLR
# but strictly speaking, the row-based centering is already done by CLR.
pca_res_A <- prcomp(clr_mat_A, center = FALSE, scale. = FALSE)

# Make a data frame of PCA scores
pca_scores_A <- as.data.frame(pca_res_A$x)
pca_scores_A$SampleID <- rownames(pca_scores_A)
# Add group info
pca_scores_A$age_group <- metadata_final$age_group

# Variation explained (PC1, PC2)
var_exp_A <- round(100 * pca_res_A$sdev^2 / sum(pca_res_A$sdev^2), 1)
xlabA <- paste0("PC1 (", var_exp_A[1], "%)")
ylabA <- paste0("PC2 (", var_exp_A[2], "%)")

# Step A5: Plot
pA <- ggplot(pca_scores_A, aes(x = PC1, y = PC2, color = age_group)) +
  geom_point(size = 2) +
  stat_ellipse(type = "norm", level = 0.75, size = 1) +
  labs(title = "PCA (CLR via compositions package)", 
       x = xlabA, 
       y = ylabA, 
       color = "Age Group") +
  theme_minimal()

# Print or save
ggsave("PCA_CLR_compositions_package.pdf", plot = pA, width = 6, height = 5)
