library(data.table)
library(DESeq2)
library(ggplot2)
library(ggfortify)

gene_count <- read.csv(file = "./PCA/gene_expression.csv")
gene_metadata <- read.csv(file = "./PCA/metadata_PCA.csv")

gene_metadata <- gene_metadata[-1,]
gene_count <- as.data.frame(t(gene_count))

names(gene_count) <- gene_count[1,]
gene_count <- gene_count[-1,]
gene_count <- sapply(gene_count, as.numeric)
gene_count <- as.data.frame(gene_count)


gene_count$condition <- gene_metadata$refinebio_disease
gene_rows <- nrow(gene_count)
gene_cols <- ncol(gene_count) -1

pca <- prcomp(gene_count[, 1:gene_cols])

pca_plot <- autoplot(pca, data = gene_count, colour = "condition")

ggsave("Gene_PCA.pdf", width = 16, height = 9)
