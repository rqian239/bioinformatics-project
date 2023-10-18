library(ggplot2)
library(ComplexHeatmap)
library(dplyr)

# Load in gene expression data
data_file <- "./data/SRP073813/SRP073813-HUGO-cleaned-top-5000.tsv"
expression_df <- readr::read_tsv(data_file) %>% tibble::column_to_rownames("first_mapped_hugo")

# Get z-scores of each row
expression_df_z <- t(apply(expression_df, 1, scale))
colnames(expression_df_z) <- colnames(expression_df)

# Load in metadata
clustering_metadata_file <- "./results/clustering/k-means-cluster-results.tsv"
clustering_metadata <- readr::read_tsv(clustering_metadata_file) %>% tibble::column_to_rownames("refinebio_accession_code")

# Group metadata by group
clustering_metadata <- clustering_metadata %>% arrange(k_cluster)

# Get order of rows when we group by k_cluster
grouping_order <- rownames(clustering_metadata)

# Reorder rows of expression_df_z
heatmap_df <- expression_df_z[, grouping_order]

# Create annotation object
anno <- HeatmapAnnotation(df = clustering_metadata, which = "column", 
    col = list(cluster = c("1" = "blue", "2" = "red", "3" = "yellow", "4" = "green", "5" = "purple"), 
    disorder_group = c("ancg_control" = "lightcoral", "ancg_schizophrenia" = "cyan")))

# Create Heatmap
Heatmap(heatmap_df, cluster_rows = T, cluster_columns = F, name = "Z-scores", bottom_annotation = anno, column_title = "Samples", row_title = "Genes")