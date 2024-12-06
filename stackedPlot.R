# stackedPlot.R
# this file contains the function
# load_data_stacked and plot_stacked_bar_plot

# Install and load necessary libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(data.table)
  library(forcats)
})

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
subset_size <- 50

# Step 2: Subset the first 50 samples and convert row names to 'Sample' column
subset_df <- sample_dataframe[1:subset_size, ] %>%
  rownames_to_column(var = "Sample")

# Step 3: Calculate total abundance for each OTU across the subset
otu_totals <- suppressWarnings(
  subset_df %>%
    select(-Sample) %>%  # Exclude the 'Sample' column
    summarise(across(everything(), sum, na.rm = TRUE)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "OTU",
      values_to = "Total_Abundance"
    )
)

# Identify the top 10 OTUs
top10_otus <- otu_totals %>%
  arrange(desc(Total_Abundance)) %>%
  slice(1:10) %>%
  pull(OTU)

# Step 4: Transform to long format
# This step is necessary for creating the stacked bar plot
# where each row represents a sample-OTU pair with the abundance value
# so we have the same sample id multiple times as each row, one for each OTU
subset_long <- subset_df %>%
  pivot_longer(
    cols = -Sample,
    names_to = "OTU",
    values_to = "Abundance"
  )
# Step 5: Categorize OTUs as Top 10 or Others
subset_long <- subset_long %>%
  mutate(OTU = ifelse(OTU %in% top10_otus, OTU, "Others"))

# Step 6: Aggregate Abundance for Each Sample and OTU
subset_long <- subset_long %>%
  group_by(Sample, OTU) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Step 7: Calculate Mean Abundance for Each OTU
otu_means <- subset_long %>%
  group_by(OTU) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Mean_Abundance))

# View OTU Means
# print(otu_means) # Uncomment to view the mean abundance of each OTU # nolint

# Step 8: Merge Mean Abundance into subset_long
subset_long <- subset_long %>%
  left_join(otu_means, by = "OTU")

# Step 9: Reorder OTU Factor Based on Mean Abundance
# This is used to reorder the OTUs in the stacked bar plot
# based on the mean abundance across all samples
# notice that the actual row will not be reordered,
# it is the factor levels that will be reordered
subset_long <- subset_long %>%
  mutate(
    OTU = fct_reorder(OTU, Mean_Abundance, .desc = TRUE),
    OTU = fct_relevel(OTU, "Others", after = Inf)
  )
# Step 10: Create the stacked bar plot
p <- ggplot(subset_long, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(
    y = "Relative Abundance",
    fill = "OTU"
  ) +
  scale_fill_brewer(palette = "Set3") +
  guides(fill = guide_legend(reverse = FALSE))

# Step 11: Save the plot
ggsave("stacked_bar_plot_top10_OTUs_ordered.pdf",
       plot = p,
       width = 14,
       height = 8)
