# UMAP Plot of Gene Expression Data

library(umap)
library(ggplot2)

# Load in data
gene_expression <- readr::read_tsv("./data/SRP073813/SRP073813-HUGO-cleaned.tsv") %>% 
  mutate(first_mapped_hugo = make.unique(as.character(first_mapped_hugo)))

metadata <- readr::read_tsv("./data/SRP073813/metadata_SRP073813_filtered.tsv") %>% tibble::column_to_rownames("refinebio_accession_code")

# Transpose the data
gene_expression_t <- as.data.frame(t(gene_expression))

# Make the first row the column names
colnames(gene_expression_t) <- gene_expression_t[1,]
gene_expression_t <- gene_expression_t[-1,]

# Make data numerics
gene_expression_t <- mutate_all(gene_expression_t, function(x) as.numeric(as.character(x)))

# Check that the metadata and gene expression data have the same rownames
all(rownames(metadata) == rownames(gene_expression_t))

umap <- umap(gene_expression_t)

# Convert UMAP results to a data frame
umap_df <- as.data.frame(umap$layout)
umap_df <- cbind(umap_df, metadata$refinebio_subject)
colnames(umap_df) <- c("UMAP_1", "UMAP_2", "disorder")

umap_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = factor(disorder))) +
  geom_point() +
  theme_minimal() +
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )


ggsave("./plots/UMAP/umap_plot.png", umap_plot, width = 16, height = 9)