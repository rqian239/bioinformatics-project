library(M3C)
library(data.table)
library(DESeq2)
library(ggplot2)
library(ggfortify)


gene_count <- read.csv(file = "./PCA/gene_expression.csv")

gene_count$hugo_id <- make.unique(gene_count$hugo_id)
gene_count <- as.data.frame(t(gene_count))
names(gene_count) <- gene_count[1,]
gene_count <- gene_count[-1,]
gene_count <- sapply(gene_count, as.numeric)
gene_count <- as.data.frame(gene_count)


gene_metadata <- read.csv(file = "./PCA/metadata_PCA.csv")
gene_metadata <- gene_metadata[-1,]

# gene_count$condition <- gene_metadata$refinebio_disease
gene_count <- unique(gene_count)

umap <- umap(gene_count)

umap_data <- data.frame(x = umap$layout[,1], y = umap$layout[,2], condition = gene_metadata$refinebio_disease)

umap_plot <- ggplot(umap_data, aes(x,y, colour = condition)) + geom_point()
ggsave("Gene_UMAP.pdf", width = 16, height = 9)
