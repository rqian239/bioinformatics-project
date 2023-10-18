library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)

gene_metadata   <- read.csv(file = "./data/SRP073813/gene_metadata.csv")
gene_metadata <- gene_metadata[-1,]

gene_expression <- read.csv(file = "./data/SRP073813/gene_expression_schizophrenia.csv")
rownames(gene_expression) <- gene_expression[,1]
gene_expression <- gene_expression[,-1]
gene_names <- colnames(gene_expression)
gene_expression <- gene_expression[order(apply(gene_expression, 1, var), decreasing = TRUE), ]
gene_expression <- gene_expression[1:5000]
gene_expression <- as.data.frame(apply(gene_expression, 1 , scale))

hclust_metadata <- read.csv(file = "./hcluster/h_cluster.csv")
hclust_metadata <- arrange(hclust_metadata, hierarchical_cluster)

# Reorder rows of expression_df_z
# Assuming you have a common identifier column in both data frames, let's call it "common_id"
common_id <- "refinebio_accession_code"  # Replace with the actual column name

# Match row names between gene_expression and hclust_metadata
gene_expression_ids <- rownames(gene_expression)
matched_rows <- match(hclust_metadata[[common_id]], gene_expression_ids)

# Remove rows with missing matches (if any)
matched_rows <- matched_rows[!is.na(matched_rows)]

# Reorder gene_expression based on the matched rows from hclust_metadata
heatmap <- t(gene_expression[matched_rows, 1:5000])
colnames(heatmap) <- gene_names[1:5000]
# Create annotation object
anno <- HeatmapAnnotation(df = hclust_metadata, which = "col", 
    col = list(hierarchical_cluster = c("1" = "blue", "2" = "red", "3" = "yellow", "4" = "green", "5" = "purple"), 
    disorder_group = c("control" = "lightcoral", "schizophrenia" = "cyan")))

# Create Heatmap
Heatmap(heatmap, cluster_rows = T, cluster_columns = F, name = "Z-scores", bottom_annotation = anno, column_title = "Samples", row_title = "Genes")
