# Script from https://alexslemonade.github.io/refinebio-examples/03-rnaseq/differential-expression_rnaseq_01.html

# Install packages if not already installed
if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}

library(DESeq2)
library(ggplot2)
library(magrittr)

# Set seed for DESeq2::plotCounts() function
set.seed(12345)

# Read in metadata TSV file
metadata_file <- "./data/SRP073813/metadata_SRP073813.tsv"
metadata <- readr::read_tsv(metadata_file)

# Read in gene expression data
data_file <- "./data/SRP073813/SRP073813-HUGO-cleaned.tsv"
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("first_mapped_hugo")
