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
                                                               "age_years",
                                                               "height_cm",
                                                               "weight_kg"))

# 5) Match sample IDs between OTU table and metadata
common_samples <- intersect(rownames(sample_dataframe_filt),
                            rownames(metadata_clean))
cat("Number of matching samples:", length(common_samples), "\n")

otu_final <- sample_dataframe_filt[common_samples, , drop = FALSE]
metadata_final <- metadata_clean[common_samples, , drop = FALSE]



# 3) Define Height groups
metadata_final <- metadata_final %>%
  mutate(height_group = ifelse(height_cm > 175, "Above175", "Below175"))

# 4) Add pseudocount and do CLR transform using compositions package
# This is a step that we need to review again later to decide for the best pseudocount value
pseudocount_value <- 1e-5
otu_version <- otu_final + pseudocount_value
otu_acomp <- acomp(otu_version)
clr_mat <- clr(otu_acomp)

# 5) PCA result (pcs_res)
pca_res <- prcomp(clr_mat, center = FALSE, scale. = FALSE)

# Extract PCA scores
pca_scores <- as.data.frame(pca_res$x)
pca_scores$SampleID <- rownames(pca_scores)

# 6) Merge PCA scores with the new `height_group`, "smoking"
pca_scores$height_group <- metadata_final[rownames(pca_scores), "height_group"]
pca_scores$smoking_frequency <- metadata_final[rownames(pca_scores), "smoking_frequency"]
pca_scores$alcohol_frequency <- metadata_final[rownames(pca_scores), "alcohol_frequency"]


# Compute % variance explained for labeling
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
xlab <- paste0("PC1 (", var_exp[1], "%)")
ylab <- paste0("PC2 (", var_exp[2], "%)")

# 7) Plot PCA with ellipses by height group
p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = height_group)) +
  geom_point(size = 0.1) +
  # stat_ellipse draws ~75% coverage ellipses if type="norm" & level=0.75
  stat_ellipse(type = "norm", level = 0.75, size = 1) +
  labs(title = "PCA (CLR) - Height Group at 175",
       x = xlab,
       y = ylab,
       color = "Height Group") +
  theme_minimal()

# 8) Save
ggsave("PCA_CLR_height175.pdf", p, width = 6, height = 5)
cat("PCA with hieght group 175 complete. Output: PCA_CLR_height175.pdf\n")

# Filter out "Not provided" entries and recode smoking responses
pca_scores_filtered <- pca_scores %>%
  filter(smoking_frequency != "Not provided") %>%
  mutate(smoking_status = ifelse(smoking_frequency == "Never", "Never", "Smoking")) %>%
  mutate(alcohol_status = ifelse(alcohol_frequency == "Never", "Never", "Alcohol"))

# Plot PCA by the new smoking_status grouping
p_smoking <- ggplot(pca_scores_filtered, aes(x = PC1, y = PC2, color = smoking_status)) +
  geom_point(size = 0.1) +
  stat_ellipse(type = "norm", level = 0.75, size = 0.1) +
  labs(title = "PCA (CLR) by Smoking Status",
       x = xlab,
       y = ylab,
       color = "Smoking Status") +
  theme_minimal()

ggsave("PCA_CLR_smoking_frequency.pdf", p_smoking, width = 6, height = 5)
cat("PCA with smoking_frequency complete. Output: PCA_CLR_smoking_frequency.pdf\n")

p_alcohol <- ggplot(pca_scores_filtered, aes(x = PC1, y = PC2, color = alcohol_status)) +
  geom_point(size = 0.1) +
  stat_ellipse(type = "norm", level = 0.75, size = 0.1) +
  labs(title = "PCA (CLR) by Alcohol Frequency",
       x = xlab,
       y = ylab,
       color = "Alcohol Frequency") +
  theme_minimal()
ggsave("PCA_CLR_alcohol_frequency.pdf", p_alcohol, width = 6, height = 5)
cat("PCA with alcohol_frequency complete. Output: PCA_CLR_alcohol_frequency.pdf\n")

pca_scores$bmi_group <- metadata_final[rownames(pca_scores), "bmi_group"]


var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
xlab <- paste0("PC1 (", var_exp[1], "%)")
ylab <- paste0("PC2 (", var_exp[2], "%)")

set.seed(123)
km_res <- kmeans(pca_scores[, c("PC1", "PC2")], centers = 3, nstart = 25)
# This vector labels each sample as cluster 1 or cluster 2
pca_scores$cluster <- factor(km_res$cluster)

p_cluster <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 0.1) +
  stat_ellipse(type = "norm", level = 0.75, size = 1) +
  labs(title = "PCA color-coded by Cluster",
       x = xlab,
       y = ylab,
       color = "Cluster") +
  theme_minimal()

ggsave("PCA_by_cluster.pdf", p_cluster, width = 6, height = 5)
cat("PCA with cluster complete. Output: PCA_by_cluster.pdf\n")

user_centers <- matrix(c(-10, 2.5,
                         -11.25, -2.5),
                       nrow = 2, byrow = TRUE)
rownames(user_centers) <- c("Center1", "Center2")
colnames(user_centers) <- c("PC1", "PC2")

user_centers

pca_scores$user_cluster <- apply(pca_scores[, c("PC1", "PC2")], 1, function(sample_coords) {
  # Distances to each user-defined center
  dists <- apply(user_centers, 1, function(center_coords) {
    sum((sample_coords - center_coords)^2)  # squared Euclidean distance
  })
  which.min(dists)  # returns the index (1 or 2) of the nearest center
})

# Convert to a factor if you prefer
pca_scores$user_cluster <- factor(pca_scores$user_cluster, levels = c(1,2))

library(ggplot2)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = user_cluster)) +
  geom_point(size = 2) +
  # (Optional) Add the centers themselves:
  geom_point(data = as.data.frame(user_centers),
             aes(x = PC1, y = PC2),
             color = "black",
             shape = 8,       # star or some distinct shape
             size = 4) +
  labs(title = "Nearest-Centroid Assignment to User-Defined Centers",
       color = "User Cluster") +
  scale_color_discrete()
theme_minimal()

# try ggpubr gg themes
# gg sci



###################################################################
# t-test
###################################################################

common_samples <- intersect(rownames(pca_scores), rownames(metadata_final))
pca_test <- pca_scores[common_samples, , drop = FALSE]
meta_test <- metadata_final[common_samples, , drop = FALSE]

# Pick which continuous variables you want to test
continuous_vars <- c("age_years", "bmi", "height_cm", "weight_kg")

for (var in continuous_vars) {
  # cluster1 vs cluster2
  values_cluster1 <- meta_test[pca_test$user_cluster == 1, var]
  values_cluster2 <- meta_test[pca_test$user_cluster == 2, var]

  # e.g. a simple t-test
  t_out <- t.test(values_cluster1, values_cluster2)

  cat("Variable:", var, "\n")
  cat("  Mean cluster1 =", mean(values_cluster1, na.rm = TRUE),
      ", Mean cluster2 =", mean(values_cluster2, na.rm = TRUE), "\n")
  cat("  t-test p-value =", t_out$p.value, "\n\n")
}



###### do the table of the relationship between pcs and metadata
# First, add the metadata variables to the PCA scores data frame.
# (Assuming these variables were cleaned in metadata_final.)
pca_scores$bmi <- metadata_final[rownames(pca_scores), "bmi"]
pca_scores$age_years <- metadata_final[rownames(pca_scores), "age_years"]
pca_scores$height_cm <- metadata_final[rownames(pca_scores), "height_cm"]
pca_scores$weight_kg <- metadata_final[rownames(pca_scores), "weight_kg"]

# Load broom package for tidying model outputs
library(broom)

# Define which principal components to analyze (e.g., PC1 to PC5)
pc_names <- grep("^PC", names(pca_scores), value = TRUE)[1:5]

# Create an empty list to store results
regression_results <- list()

# Loop through each PC, fit a linear model with metadata variables, and store results
for (pc in pc_names) {
  # adjust the formula later
  formula <- as.formula(paste(pc, "~ bmi + age_years + height_cm + weight_kg + smoking_frequency + alcohol_frequency"))

  # Fit the linear model
  fit <- lm(formula, data = pca_scores)

  # Tidy up the results to extract beta and p-values for each predictor.
  tidy_fit <- tidy(fit)

  # Also extract the overall R-squared from the model summary.
  r_squared <- summary(fit)$r.squared

  # Add a column for the PC name and overall R-squared.
  tidy_fit$PC <- pc
  tidy_fit$R_squared <- r_squared

  # Store the results
  regression_results[[pc]] <- tidy_fit
}

# Combine all results into a single data frame
regression_df <- do.call(rbind, regression_results)

# View the results
print(regression_df)
