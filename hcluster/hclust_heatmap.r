library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)

top_variable_genes <- 500

gene_metadata <- read.csv(file = "./data/SRP073813/gene_metadata.csv")
gene_metadata <- gene_metadata[-1, ]

gene_expression <- read.csv(file = "./data/SRP073813/gene_expression_schizophrenia.csv")
gene_names <- colnames(gene_expression)
row.names(gene_expression) <- gene_expression[, 1]
gene_expression <- gene_expression[, -1]
gene_expression <- as.data.frame(apply(gene_expression, 1, scale))

gene_names <- gene_names[-1]
row.names(gene_expression) <- gene_names

gene_expression <- gene_expression[order(apply(gene_expression, 1, var), decreasing = TRUE), ]
gene_expression <- gene_expression[1:top_variable_genes, ]
gene_expression <- t(gene_expression)

hclust_metadata <- read.csv(file = "./hcluster/h_cluster.csv")
hclust_metadata <- arrange(hclust_metadata, hierarchical_cluster)

# Reorder rows of expression_df_z
# Assuming you have a common identifier column in both data frames, let's call it "common_id"
common_id <- "refinebio_accession_code" # Replace with the actual column name

# Match row names between gene_expression and hclust_metadata
gene_expression_ids <- rownames(gene_expression)
matched_rows <- match(hclust_metadata[[common_id]], gene_expression_ids)

# Reorder gene_expression based on the matched rows from hclust_metadata
heatmap <- gene_expression[matched_rows, ]
# Create annotation object
anno <- HeatmapAnnotation(
    df = hclust_metadata, which = "col",
    col = list(
        hierarchical_cluster = c("1" = "blue", "2" = "red", "3" = "yellow", "4" = "green", "5" = "purple"),
        disorder_group = c("control" = "lightcoral", "schizophrenia" = "cyan")
    )
)


# Create Heatmap
# Heatmap(heatmap, cluster_rows = T, cluster_columns = F, name = "Z-scores", bottom_annotation = anno, column_title = "Samples", row_title = "Genes")
filename <- sprintf("./hcluster/Hierarchical_Heatmap_%d.pdf", top_variable_genes)
pdf(filename, width = 9, height = 16)
heatmap_plot <- Heatmap(t(heatmap), cluster_rows = T, cluster_columns = F, name = "Z-scores", bottom_annotation = anno, column_title = "Samples", row_title = "Genes")

dev.off()
