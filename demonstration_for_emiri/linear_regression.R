# load libraries

library(data.table)
library(CCA)
library(CCP)
library(dplyr)
library(ggplot2)

# load the data


otu_file <- "subset_tables/subset_table_order.tsv"
source("load_data.R")
sample_dataframe <- load_data_stacked(otu_file)

cutoff <- 0.01
sample_dataframe_filt <- filter_low_abundance_taxa(sample_dataframe, cutoff)

metadata_file <- "dataset/metadata.tsv"
metadata_clean <- load_and_filter_metadata(metadata_file,
                                           sample_col_name = "sample_name",
                                           continuous_cols = c("bmi",
                                                               "age_years"))
common_samples <- intersect(rownames(sample_dataframe_filt),
                            rownames(metadata_clean))
cat("Number of matching samples:", length(common_samples), "\n")

otu_final <- sample_dataframe_filt[common_samples, , drop = FALSE]
metadata_final <- metadata_clean[common_samples, , drop = FALSE]

bifido_abundance <- otu_final[,
                              "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Bifidobacteriales"] #nolint



df_model <- data.frame(
  age_years        = metadata_final$age_years,
  bifido_abundance = bifido_abundance
)

model_simple <- lm(bifido_abundance ~ age_years, data = df_model)
summary(model_simple)
p <- ggplot(df_model, aes(x = age_years, y = bifido_abundance)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Bifidobacterium abundance vs. age_years",
       x = "Age (years)",
       y = "Bifidobacterium abundance")
    
print(p)
