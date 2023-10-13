# Perform k-means clustering on gene expression data
# From https://www.youtube.com/watch?v=NKQpVU1LTm8

# install.packages("factoextra")

# Load libraries
library(factoextra)
library(dplyr)

# Read in gene expression data
gene_expression <- readr::read_tsv("./data/SRP073813/SRP073813-HUGO-cleaned-top-5000.tsv")

# Transpose the data
gene_expression_t <- t(gene_expression)

# Make the first row the column names
colnames(gene_expression_t) <- gene_expression_t[1,]
gene_expression_t <- gene_expression_t[-1,]

# Convert matrix to data frame
gene_expression_t <- as.data.frame(gene_expression_t)


# Load in metadata
metadata <- readr::read_tsv("./data/SRP073813/metadata_SRP073813_filtered.tsv") %>% tibble::column_to_rownames("refinebio_accession_code")
# Check rownames match
all(rownames(metadata) == rownames(gene_expression_t))

# Get labels of samples
sample_labels <- metadata$refinebio_subject
table(sample_labels)


# Scale the data
gene_expression_t_scaled <- scale(gene_expression_t)