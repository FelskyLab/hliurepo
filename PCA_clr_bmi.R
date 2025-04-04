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


# Extract loadings for each feature (OTU)
pca_res_A <- pca_res
pca_loadings <- as.data.frame(pca_res_A$rotation)
pca_loadings$Feature <- rownames(pca_loadings)  # keep the OTU names in a column

# Variation explained, for axis labels
var_exp_A <- round(100 * pca_res_A$sdev^2 / sum(pca_res_A$sdev^2), 1)
xlabA <- paste0("PC1 (", var_exp_A[1], "%)")
ylabA <- paste0("PC2 (", var_exp_A[2], "%)")

# Decide on a scale factor for the loadings arrows
# A common quick approach is to pick a factor so that arrows
# fit roughly within the sample cloud. You can tweak or compute it dynamically.
scaleFactor <- 5

# Create new columns for the scaled loadings
pca_loadings$PC1_arrow <- pca_loadings$PC1 * scaleFactor
pca_loadings$PC2_arrow <- pca_loadings$PC2 * scaleFactor

pca_scores_A <- pca_scores
# Base plot using sample scores
pBiplot <- ggplot(pca_scores_A, aes(x = PC1, y = PC2, color = bmi_group)) +
  geom_point(size = 2) +
  stat_ellipse(type = "norm", level = 0.75, size = 1) +
  labs(title = "PCA Biplot (CLR) with OTU Loadings",
       x = xlabA,
       y = ylabA,
       color = "Age Group") +
  theme_minimal()

# Add loadings as arrows (geom_segment) + label them with geom_text
pBiplot <- pBiplot +
  geom_segment(data = pca_loadings,
               aes(x = 0, y = 0,
                   xend = PC1_arrow, yend = PC2_arrow),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black",
               alpha = 0.7) +
  geom_text(data = pca_loadings,
            aes(x = PC1_arrow, y = PC2_arrow, label = Feature),
            color = "black",
            size = 3,
            nudge_y = 0.0)

# Print or save
print(pBiplot)
ggsave("PCA_Loadings_Biplot.pdf", pBiplot, width = 6, height = 5)