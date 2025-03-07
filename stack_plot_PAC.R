# stackedPlot.R
# This file shows how to keep only the order (o__) and remove everything else.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(data.table)
  library(forcats)
})

# --------------------------------------------------------------------
# Step 1: Load data function
# --------------------------------------------------------------------
load_data_stacked <- function(file_path) {
  mydata <- fread(file_path, header = TRUE, sep = "\t", skip = 1)
  cat("Original data dimensions: ", dim(mydata), "\n")

  otu_ids <- mydata[[1]]
  otu_counts <- as.matrix(mydata[, -1])
  rownames(otu_counts) <- otu_ids

  # Transpose so that samples are rows
  sample_data <- t(otu_counts)
  sample_dataframe <- as.data.frame(sample_data)
  return(sample_dataframe)
}

# --------------------------------------------------------------------
# Step 2: Load and subset data
# --------------------------------------------------------------------
sample_dataframe <- load_data_stacked("subset_tables/subset_table_order.tsv")
subset_size <- 1000

subset_df <- sample_dataframe[1:subset_size, ] %>%
  rownames_to_column(var = "Sample")

# --------------------------------------------------------------------
# Step 3: Identify top 10 OTUs
# --------------------------------------------------------------------
otu_totals <- suppressWarnings(
  subset_df %>%
    select(-Sample) %>%
    summarise(across(everything(), sum, na.rm = TRUE)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "OTU",
      values_to = "Total_Abundance"
    )
)

top10_otus <- otu_totals %>%
  arrange(desc(Total_Abundance)) %>%
  slice(1:10) %>%
  pull(OTU)

# --------------------------------------------------------------------
# Step 4: Convert to long format
# --------------------------------------------------------------------
subset_long <- subset_df %>%
  pivot_longer(
    cols = -Sample,
    names_to = "OTU",
    values_to = "Abundance"
  )

# --------------------------------------------------------------------
# Step 5: Categorize OTUs as Top 10 or Others
# --------------------------------------------------------------------
subset_long <- subset_long %>%
  mutate(OTU = ifelse(OTU %in% top10_otus, OTU, "Others"))

# --------------------------------------------------------------------
# Step 6: Aggregate abundance for each sample and OTU
# --------------------------------------------------------------------
subset_long <- subset_long %>%
  group_by(Sample, OTU) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# --------------------------------------------------------------------
# Step 7: Calculate mean abundance per OTU
# --------------------------------------------------------------------
otu_means <- subset_long %>%
  group_by(OTU) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Mean_Abundance))

# --------------------------------------------------------------------
# Step 8: Merge mean abundance into subset_long
# --------------------------------------------------------------------
subset_long <- subset_long %>%
  left_join(otu_means, by = "OTU")

# --------------------------------------------------------------------
# Step 9: Reorder OTU factor based on mean abundance
# --------------------------------------------------------------------
subset_long <- subset_long %>%
  mutate(
    OTU = fct_reorder(OTU, Mean_Abundance, .desc = TRUE),
    OTU = fct_relevel(OTU, "Others", after = Inf)
  )

# --------------------------------------------------------------------
# Step 9b: Keep only the order portion (o__)
# --------------------------------------------------------------------
# We extract the text starting at 'o__' and continuing until the next semicolon
# (or the string's end), ignoring everything else.
# For 'Others', we keep it as 'Others'.
subset_long <- subset_long %>%
  mutate(OTU = ifelse(
    OTU != "Others",
    sub(".*(o__[^;]+).*", "\\1", OTU),  # keep only 'o__...'
    "Others"
  ))

# --------------------------------------------------------------------
# Step 10: Reorder samples based on top OTU abundance
# --------------------------------------------------------------------
top_otu <- top10_otus[1]

# If you want to also convert the top_otu itself to its shortened form:
top_otu_short <- sub(".*(o__[^;]+).*", "\\1", top_otu)

top_otu_abundance <- subset_long %>%
  filter(OTU == top_otu_short) %>%
  select(Sample, Abundance)

sample_order <- top_otu_abundance %>%
  arrange(desc(Abundance)) %>%
  pull(Sample)

subset_long <- subset_long %>%
  mutate(Sample = factor(Sample, levels = sample_order))

# --------------------------------------------------------------------
# Step 11: Create the stacked bar plot
# --------------------------------------------------------------------
p <- ggplot(subset_long, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", width = 0.2) +
  theme_bw() +
  theme(
    # Remove sample names on the x-axis
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),

    # Increase axis text sizes
    axis.text.y  = element_text(size = 14),
    axis.title.y = element_text(size = 16),

    # Optionally increase legend font size
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  labs(
    y = "Relative Abundance",
    fill = "Order"
  ) +
  scale_fill_brewer(palette = "Set3") +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_x_discrete(limits = sample_order) +
  coord_cartesian(clip = "off")

ggsave("stacked_bar_plot_top10_OTUs_ordered_new.pdf",
       plot = p,
       width = 40,
       height = 8)