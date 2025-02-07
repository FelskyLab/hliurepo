#!/usr/bin/env Rscript

#  Example R Script: Loading, Filtering, Excluding Missing Y Data, 
#                    and Running CCA
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

#################################################################
# 3) Load the OTU table and filter low-abundance taxa
#################################################################

# Path to your OTU table (example)
otu_file <- "subset_tables/subset_table_order.tsv"

# Load and transpose
sample_dataframe <- load_data_stacked(otu_file)
cat("After transpose, dimension is: ", dim(sample_dataframe), "\n")
# Now 'sample_dataframe' has rows = samples, columns = taxa/features

# Filter out taxa with average abundance < 0.1% (0.001 in fractional form)
threshold <- 0.01  # 1%
col_means <- colMeans(sample_dataframe)
sample_dataframe_filt <- sample_dataframe[, col_means >= threshold, drop = FALSE]

cat("Dimension after filtering (< 0.1% avg): ", dim(sample_dataframe_filt), "\n")

#################################################################
# 4) Load metadata and exclude samples with missing BMI/age
#################################################################

# Suppose you have a metadata file with sample_name, BMI, age_years, etc.
metadata_file <- "dataset/metadata.tsv"
metadata <- fread(metadata_file, header = TRUE, sep = "\t") %>%
  as.data.frame()

# For example, you might have:
#   sample_name   bmi    age_years
#   Sample1       22.5   34
#   Sample2       27.1   56
#   ...

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


### 9000 samples overlap between OTU and metadata ###

# Subset both data frames to matching samples
otu_final <- sample_dataframe_filt[common_samples, , drop=FALSE]
metadata_final <- metadata_clean[common_samples, , drop=FALSE]


# 5) Prepare blocks for Canonical Correlation Analysis (CCA)
# Suppose we want:
#    X = microbiome features (columns in otu_final)
#    Y = BMI + age_years (columns in metadata_final)

X <- otu_final
Y <- metadata_final[, c("bmi", "age_years")]

common_samples <- intersect(rownames(X), rownames(Y))
common_samples_sorted <- sort(common_samples)

X_sub <- X[common_samples_sorted, , drop = FALSE]
Y_sub <- Y[common_samples_sorted, , drop = FALSE]

# Optional: Scale (mean=0, sd=1)
X_scaled <- scale(X_sub, center = TRUE, scale = TRUE)
Y_scaled <- scale(Y_sub, center = TRUE, scale = TRUE)


# 6) Run CCA

res.cc <- cc(X_scaled, Y_scaled)

# Inspect the canonical correlations
cat("Canonical correlations:\n")
print(res.cc$cor)

# Canonical coefficients (loadings) for X (bacteria side)
cat("X side coefficients (first few):\n")
print(head(res.cc$xcoef))

# Canonical coefficients (loadings) for Y (metadata side)
cat("Y side coefficients:\n")
print(res.cc$ycoef)

########################################################################
# 7) (Optional) Significance Testing
########################################################################

# Using the CCP package to get approximate p-values
p_values <- p.asym(res.cc$cor, nrow(X_scaled), ncol(X_scaled), ncol(Y_scaled))
cat("CCA dimension significance:\n")
print(p_values)

########################################################################
# 8) Next Steps
########################################################################

# - Interpret which taxa have the largest loadings in the first canonical dimension(s).
# - Check which dimension is significant (p-value).
# - Possibly refine the threshold for filtering or add more metadata variables.
# - Explore partial correlation or partial least squares if needed.
# - Visualize canonical variates if desired.

# Extract the taxon abundance for Bifidobacteriales
bifido_abundance <- otu_final[ , "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Bifidobacteriales"]
df_model <- data.frame(
  age_years        = metadata_final$age_years,
  bmi              = metadata_final$bmi,
  bifido_abundance = bifido_abundance
)

head(df_model)

model_simple <- lm(bifido_abundance ~ age_years, data = df_model)
summary(model_simple)

# 1) Create a basic scatter plot
library(ggplot2)

p1 <- ggplot(df_model, aes(x = age_years, y = bifido_abundance)) +
        geom_point(
          shape = 16,  # small solid circle
          size = 1,    # smaller point size
          alpha = 0.8, # transparency for overlapping points
          color = "blue"
        ) +
        geom_smooth(method = "lm", color = "red", se = TRUE) +
        labs(
          title = "Scatter Plot of Bifidobacteriales Abundance vs. Age",
          x = "Age (years)",
          y = "Bifidobacteriales Abundance"
        ) +
        theme_minimal()

ggsave("scatter_plot_bifido_age.pdf", plot = p1, width = 6, height = 4, dpi = 300)

# Create a new column that is log-transformed
df_model$log_bifido_abundance <- log1p(df_model$bifido_abundance)

# Now fit the linear model using the transformed abundance
model_simple_log <- lm(log_bifido_abundance ~ age_years, data = df_model)
summary(model_simple_log)

p2 <- ggplot(df_model, aes(x = age_years, y = bifido_abundance)) +
        geom_point(size = 1, alpha = 0.8, shape = 16, color = "blue") +
        geom_smooth(method = "lm", color = "red", se = TRUE) +
        scale_y_log10() +
        labs(
          title = "Bifido Abundance vs. Age (log-scale y-axis)",
          x = "Age (years)",
          y = "Bifidobacteriales Abundance (log scale)"
        ) +
        theme_minimal()

ggsave("scatter_plot_bifido_age_log.pdf", plot = p2, width = 6, height = 4, dpi = 300)