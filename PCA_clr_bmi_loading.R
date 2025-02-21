# 1) Load libraries
library(dplyr)
library(ggplot2)
library(compositions)

# 2) Source your helper scripts
source("load_data.R")        # Contains load_data_stacked()

# 3) Load and filter the OTU table
otu_file <- "subset_tables/subset_table_order.tsv"
sample_dataframe <- load_data_stacked(otu_file)

# e.g. filter out taxa with < 1% average abundance
threshold <- 0.01
sample_dataframe_filt <- filter_low_abundance_taxa(sample_dataframe, threshold)

# 4) Load and filter the metadata
metadata_file <- "dataset/metadata.tsv"
# e.g., we want 'bmi' and 'age_years' to be numeric and drop NA
metadata_clean <- load_and_filter_metadata(metadata_file,
                                           sample_col_name = "sample_name",
                                           continuous_cols = c("bmi",
                                                               "age_years"))

# 5) Match sample IDs between OTU table and metadata
common_samples <- intersect(rownames(sample_dataframe_filt),
                            rownames(metadata_clean))
cat("Number of matching samples:", length(common_samples), "\n")

otu_final <- sample_dataframe_filt[common_samples, , drop = FALSE]
metadata_final <- metadata_clean[common_samples, , drop = FALSE]



# 3) Define BMI groups
metadata_final <- metadata_final %>%
  mutate(bmi_group = ifelse(bmi > 24, "Above24", "Below24"))

# 4) (Optional) Add pseudocount and do CLR transform using compositions package
pseudocount_value <- 1e-5
otu_version <- otu_final + pseudocount_value
otu_acomp <- acomp(otu_version)
clr_mat <- clr(otu_acomp)

# 5) PCA
pca_res <- prcomp(clr_mat, center = FALSE, scale. = FALSE)

# Extract PCA scores
pca_scores <- as.data.frame(pca_res$x)
pca_scores$SampleID <- rownames(pca_scores)

# 6) Merge PCA scores with the new `bmi_group`
pca_scores$bmi_group <- metadata_final[rownames(pca_scores), "bmi_group"]

# Compute % variance explained for labeling
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
xlab <- paste0("PC1 (", var_exp[1], "%)")
ylab <- paste0("PC2 (", var_exp[2], "%)")

# 7) Plot PCA with ellipses by BMI group
p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = bmi_group)) +
  geom_point(size = 0.1) +
  # stat_ellipse draws ~75% coverage ellipses if type="norm" & level=0.75
  stat_ellipse(type = "norm", level = 0.75, size = 1) +
  labs(title = "PCA (CLR) - BMI grouping at 24",
       x = xlab,
       y = ylab,
       color = "BMI Group") +
  theme_minimal()

# 8) Save
ggsave("PCA_CLR_bmi24.pdf", p, width = 6, height = 5)
cat("PCA with BMI grouping (24) complete. Output: PCA_CLR_bmi24.pdf\n")

# Extract loadings for each OTU
pca_loadings <- as.data.frame(pca_res$rotation)
pca_loadings$Taxonomy <- rownames(pca_loadings)  # keep the OTU IDs/labels

# Example only if you have a separate mapping
# pca_loadings <- left_join(pca_loadings, taxonomy_map, by = "OTU")

scaleFactor <- 5
pca_loadings$PC1_scaled <- pca_loadings$PC1 * scaleFactor
pca_loadings$PC2_scaled <- pca_loadings$PC2 * scaleFactor

# Variation explained (for axis labels)
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
xlab <- paste0("PC1 (", var_exp[1], "%)")
ylab <- paste0("PC2 (", var_exp[2], "%)")

p_tax_legend <- ggplot(pca_loadings, aes(x = PC1_scaled, y = PC2_scaled, color = Taxonomy)) +
  geom_point(size = 3) +
  labs(title = "PCA Loadings by Full Taxonomy (Legend)",
       x = xlab,
       y = ylab,
       color = "Full Taxonomy") +
  theme_minimal()

ggsave("PCA_Loadings_FullTax_Legend.pdf", p_tax_legend, width = 8, height = 6)
