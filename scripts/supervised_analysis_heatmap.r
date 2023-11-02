library(ggplot2)
library(ComplexHeatmap)
library(dplyr)

num_genes_in_dataset <- 5000

# Load in gene expression data
data_file <- sprintf("./data/SRP073813/SRP073813-HUGO-cleaned-top-%d.tsv", num_genes_in_dataset)
expression_df <- readr::read_tsv(data_file) %>% tibble::column_to_rownames("first_mapped_hugo")

# Load in gene signature data
gene_signatures_filepath <- sprintf("./results/svm/linear-svm-%d-gene-signatures.tsv", num_genes_in_dataset)
gene_signatures_df <- readr::read_tsv(gene_signatures_filepath)

# Only include genes that are in the gene signature
expression_df <- expression_df[rownames(expression_df) %in% gene_signatures_df$gene_signatures, ]



# Load in metadata
total_metadata <- readr::read_tsv("./data/SRP073813/metadata_SRP073813_filtered.tsv") %>% tibble::column_to_rownames("refinebio_accession_code")
predictions_filepath <- sprintf("./results/svm/linear-svm-%d-predictions.tsv", num_genes_in_dataset)
predictions_df <- readr::read_tsv(predictions_filepath) %>% tibble::column_to_rownames("Gene")

# Check rownames match
all(rownames(total_metadata) == rownames(predictions_df))

# Combine metadata and predictions
metadata_for_annotations <- cbind(disorder = total_metadata$refinebio_subject, predictions_df)

# Get z-scores of each row
expression_df_z <- t(apply(expression_df, 1, scale))
colnames(expression_df_z) <- colnames(expression_df)


# Group metadata by group
metadata_for_annotations <- metadata_for_annotations %>% arrange(disorder)

# Get order of rows when we group by k_cluster
grouping_order <- rownames(metadata_for_annotations)

# Reorder rows of expression_df_z
heatmap_df <- expression_df_z[, grouping_order]

# Create annotation object
anno <- HeatmapAnnotation(df = metadata_for_annotations, which = "column", 
    col = list(disorder = c("ancg_control" = "lightcoral", "ancg_schizophrenia" = "cyan"), 
    predictions = c("ancg_control" = "red", "ancg_schizophrenia" = "blue")))

# Create Heatmap
Heatmap(heatmap_df, cluster_rows = T, cluster_columns = F, name = "Z-scores", bottom_annotation = anno, column_title = "Samples", row_title = "Genes")