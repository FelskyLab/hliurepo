# Description: Load OTU table data from a file and return as a data frame.
# rownames are samples, columns are taxa.
library(data.table)

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

library(dplyr)

# 2.1) Filter out taxa/features that have low average abundance
#      threshold is assumed to be a fraction (e.g., 0.01 for 1%).
filter_low_abundance_taxa <- function(sample_df, threshold = 0.01) {
  col_means <- colMeans(sample_df)
  # Keep columns whose mean >= threshold
  filtered_df <- sample_df[, col_means >= threshold, drop = FALSE]
  cat("Dimension after filtering (<", threshold, "avg):", dim(filtered_df), "\n")
  return(filtered_df)
}


# 2.2) Load metadata, filter out invalid numeric columns for continuous variables
load_and_filter_metadata <- function(metadata_file,
                                     sample_col_name = "sample_name",
                                     continuous_cols = c("bmi",
                                                         "age_years",
                                                         "height_cm",
                                                         "weight_kg")) {
  meta <- fread(metadata_file, header = TRUE, sep = "\t") %>%
    as.data.frame()

  # Make sure sample_col_name is present
  if (!sample_col_name %in% names(meta)) {
    stop("Metadata is missing the specified sample column: ", sample_col_name)
  }

  # Set rownames to the sample identifier
  rownames(meta) <- meta[[sample_col_name]]

  # Convert columns of interest to numeric, coerce invalid strings to NA
  for (col in continuous_cols) {
    if (col %in% colnames(meta)) {
      meta[[col]] <- as.numeric(as.character(meta[[col]]))
    } else {
      warning("Column ", col, " not found in metadata; ignoring.")
    }
  }

  # Filter out rows with NA in any of the requested continuous columns
  meta_clean <- meta %>%
    filter_at(vars(continuous_cols), all_vars(!is.na(.)))

  cat("Metadata dimension after removing NA in", paste(continuous_cols, collapse=", "), #nolint
      ":", dim(meta_clean), "\n")

  return(meta_clean)
}