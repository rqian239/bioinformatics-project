# This script will sort the top X genes by variance

library(dplyr)
library(tibble)

X <- 5000

# Load in gene expression data
data_file <- "./data/SRP073813/SRP073813-HUGO-cleaned.tsv"
expression_df <- readr::read_tsv(data_file) %>% 
  mutate(first_mapped_hugo = make.unique(as.character(first_mapped_hugo))) %>%
  tibble::column_to_rownames("first_mapped_hugo")

# Sort by variance
sorted_expression_df <- expression_df[order(apply(expression_df, 1, var), decreasing = TRUE), ]

# # CHECK: Add a variance row to check that the sorting worked
# sorted_expression_df <- sorted_expression_df %>% 
#   mutate(variance = apply(sorted_expression_df, 1, var))

# Include only the top X genes by variance
top_X_genes <- sorted_expression_df[1:X, ]

# Convert rownames back to a column
output_df <- tibble::rownames_to_column(top_X_genes, "first_mapped_hugo")

# Write to a TSV
tsv_filename <- sprintf("SRP073813-HUGO-cleaned-top-%d.tsv", X)
write.table(output_df, tsv_filename, sep = "\t", quote = FALSE, row.names = FALSE)