# 1) Load libraries
library(dplyr)
library(ggplot2)
library(compositions)
library(broom)  # for tidying model outputs

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
# e.g., we want 'bmi', 'age_years', 'height_cm', and 'weight_kg' to be numeric and drop NA
metadata_clean <- load_and_filter_metadata(metadata_file,
                                           sample_col_name = "sample_name",
                                           continuous_cols = c("bmi",
                                                               "age_years",
                                                               "height_cm",
                                                               "weight_kg"))

# 5) Remove extreme values from metadata
metadata_clean <- metadata_clean %>%
  filter(bmi <= 50,
         height_cm >= 100,
         weight_kg <= 200,
         age_years <= 110)

# 6) Match sample IDs between OTU table and metadata
common_samples <- intersect(rownames(sample_dataframe_filt), rownames(metadata_clean))
cat("Number of matching samples:", length(common_samples), "\n")
otu_final <- sample_dataframe_filt[common_samples, , drop = FALSE]
metadata_final <- metadata_clean[common_samples, , drop = FALSE]

# 7) Define Height groups in metadata
metadata_final <- metadata_final %>%
  mutate(height_group = ifelse(height_cm > 175, "Above175", "Below175"))

# 8) Add pseudocount and do CLR transform using the compositions package
# (This step may be reviewed later to decide on the best pseudocount value)
pseudocount_value <- 1e-5
otu_version <- otu_final + pseudocount_value
otu_acomp <- acomp(otu_version)
clr_mat <- clr(otu_acomp)

# 9) Compute PCA from the CLR-transformed data
pca_res <- prcomp(clr_mat, center = FALSE, scale. = FALSE)
# Extract PCA scores and add sample IDs
pca_scores <- as.data.frame(pca_res$x)
pca_scores$SampleID <- rownames(pca_scores)

# 10) Merge PCA scores with metadata variables for plotting
pca_scores$height_group <- metadata_final[rownames(pca_scores), "height_group"]
pca_scores$smoking_frequency <- metadata_final[rownames(pca_scores), "smoking_frequency"]
pca_scores$alcohol_frequency <- metadata_final[rownames(pca_scores), "alcohol_frequency"]

# 11) Compute % variance explained for labeling the axes
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
xlab <- paste0("PC1 (", var_exp[1], "%)")
ylab <- paste0("PC2 (", var_exp[2], "%)")

# 12) Plot PCA with ellipses by height group and save
p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = height_group)) +
  geom_point(size = 0.1) +
  stat_ellipse(type = "norm", level = 0.75, size = 1) +
  labs(title = "PCA (CLR) - Height Group at 175",
       x = xlab,
       y = ylab,
       color = "Height Group") +
  theme_minimal()
ggsave("PCA_CLR_height175.pdf", p, width = 6, height = 5)
cat("PCA with height group 175 complete. Output: PCA_CLR_height175.pdf\n")

# 13) Filter out "Not provided" entries and recode smoking and alcohol responses
pca_scores_filtered <- pca_scores %>%
  filter(smoking_frequency != "Not provided") %>%
  mutate(smoking_status = ifelse(smoking_frequency == "Never", "Never", "Smoking"),
         alcohol_status = ifelse(alcohol_frequency == "Never", "Never", "Alcohol"))

# 14) Plot PCA by smoking status and save
p_smoking <- ggplot(pca_scores_filtered, aes(x = PC1, y = PC2, color = smoking_status)) +
  geom_point(size = 0.1) +
  stat_ellipse(type = "norm", level = 0.75, size = 0.1) +
  labs(title = "PCA (CLR) by Smoking Status",
       x = xlab,
       y = ylab,
       color = "Smoking Status") +
  theme_minimal()
ggsave("PCA_CLR_smoking_frequency.pdf", p_smoking, width = 6, height = 5)
cat("PCA with smoking status complete. Output: PCA_CLR_smoking_frequency.pdf\n")

# 15) Plot PCA by alcohol status and save
p_alcohol <- ggplot(pca_scores_filtered, aes(x = PC1, y = PC2, color = alcohol_status)) +
  geom_point(size = 0.1) +
  stat_ellipse(type = "norm", level = 0.75, size = 0.1) +
  labs(title = "PCA (CLR) by Alcohol Frequency",
       x = xlab,
       y = ylab,
       color = "Alcohol Frequency") +
  theme_minimal()
ggsave("PCA_CLR_alcohol_frequency.pdf", p_alcohol, width = 6, height = 5)
cat("PCA with alcohol frequency complete. Output: PCA_CLR_alcohol_frequency.pdf\n")

# 16) Assign BMI group (if available in metadata)
pca_scores$bmi_group <- metadata_final[rownames(pca_scores), "bmi_group"]

# 17) Recompute axis labels (if needed)
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
xlab <- paste0("PC1 (", var_exp[1], "%)")
ylab <- paste0("PC2 (", var_exp[2], "%)")

# 18) K-means clustering on the first two PCs
set.seed(123)
km_res <- kmeans(pca_scores[, c("PC1", "PC2")], centers = 3, nstart = 25)
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

# 19) Nearest-centroid assignment to user-defined centers
user_centers <- matrix(c(-10, 2.5,
                         -11.25, -2.5),
                       nrow = 2, byrow = TRUE)
rownames(user_centers) <- c("Center1", "Center2")
colnames(user_centers) <- c("PC1", "PC2")
print(user_centers)

pca_scores$user_cluster <- apply(pca_scores[, c("PC1", "PC2")], 1, function(sample_coords) {
  # Compute squared Euclidean distances to each user-defined center
  dists <- apply(user_centers, 1, function(center_coords) {
    sum((sample_coords - center_coords)^2)
  })
  which.min(dists)  # returns the index (1 or 2) of the nearest center
})
pca_scores$user_cluster <- factor(pca_scores$user_cluster, levels = c(1,2))

ggplot(pca_scores, aes(x = PC1, y = PC2, color = user_cluster)) +
  geom_point(size = 2) +
  # (Optional) Add the centers themselves:
  geom_point(data = as.data.frame(user_centers),
             aes(x = PC1, y = PC2),
             color = "black",
             shape = 8,       # distinct shape for centers
             size = 4) +
  labs(title = "Nearest-Centroid Assignment to User-Defined Centers",
       color = "User Cluster") +
  scale_color_discrete() +
  theme_minimal()

# 20) t-test comparing continuous variables between user clusters
common_samples <- intersect(rownames(pca_scores), rownames(metadata_final))
pca_test <- pca_scores[common_samples, , drop = FALSE]
meta_test <- metadata_final[common_samples, , drop = FALSE]
continuous_vars <- c("age_years", "bmi", "height_cm", "weight_kg")

for (var in continuous_vars) {
  # cluster1 vs cluster2
  values_cluster1 <- meta_test[pca_test$user_cluster == 1, var]
  values_cluster2 <- meta_test[pca_test$user_cluster == 2, var]

  t_out <- t.test(values_cluster1, values_cluster2)

  cat("Variable:", var, "\n")
  cat("  Mean cluster1 =", mean(values_cluster1, na.rm = TRUE),
      ", Mean cluster2 =", mean(values_cluster2, na.rm = TRUE), "\n")
  cat("  t-test p-value =", t_out$p.value, "\n\n")
}

# 21) Analyze the relationship between principal components and metadata
# Add continuous metadata variables to the PCA scores data frame
pca_scores$bmi <- metadata_final[rownames(pca_scores), "bmi"]
pca_scores$age_years <- metadata_final[rownames(pca_scores), "age_years"]
pca_scores$height_cm <- metadata_final[rownames(pca_scores), "height_cm"]
pca_scores$weight_kg <- metadata_final[rownames(pca_scores), "weight_kg"]

# Define which principal components to analyze (e.g., PC1 to PC5)
pc_names <- grep("^PC", names(pca_scores), value = TRUE)[1:5]
regression_results <- list()

# Loop through each PC, fit a linear model with metadata predictors, and store results
for (pc in pc_names) {
  formula <- as.formula(paste(pc, "~ bmi + age_years + height_cm + weight_kg + smoking_frequency + alcohol_frequency"))
  fit <- lm(formula, data = pca_scores)
  tidy_fit <- tidy(fit)
  r_squared <- summary(fit)$r.squared
  tidy_fit$PC <- pc
  tidy_fit$R_squared <- r_squared
  regression_results[[pc]] <- tidy_fit
}

# Combine regression results into a single data frame and print
regression_df <- do.call(rbind, regression_results)
print(regression_df)


# create a elbow plot
# Create a data frame for the variance explained by each PC
num_pcs <- length(pca_res$sdev)
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pc_df <- data.frame(PC = 1:num_pcs, VarianceExplained = var_exp)

#22) ### Create an elbow plot using ggplot2
p_elbow <- ggplot(pc_df, aes(x = PC, y = VarianceExplained)) +
  geom_line(group = 1) +
  geom_point(size = 2) +
  labs(title = "PCA Elbow Plot",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal()

# Save the elbow plot
ggsave("PCA_Elbow_Plot.pdf", p_elbow, width = 6, height = 5)
cat("PCA Elbow Plot complete. Output: PCA_Elbow_Plot.pdf\n")

# 23) Geography separation: Plot PCA by country of birth
# Merge country_of_birth from metadata into pca_scores
pca_scores$country_of_birth <- metadata_final[rownames(pca_scores), "country_of_birth"]

# Filter for the desired countries
pca_scores_geo <- pca_scores %>% 
  filter(country_of_birth %in% c("United States", "United Kingdom", "Australia"))

p_geo <- ggplot(pca_scores_geo, aes(x = PC1, y = PC2, color = country_of_birth)) +
  geom_point(size = 0.1) +
  stat_ellipse(type = "norm", level = 0.75, size = 1) +
  labs(title = "PCA (CLR) by Country of Birth",
       x = xlab,
       y = ylab,
       color = "Country of Birth") +
  theme_minimal()

ggsave("PCA_CLR_country_of_birth.pdf", p_geo, width = 6, height = 5)
cat("PCA with country of birth separation complete. Output: PCA_CLR_country_of_birth.pdf\n")