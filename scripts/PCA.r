# PCA plot on gene expression data
# Code adapted from Dr. Bharatendra Rai: https://www.kaggle.com/code/bharatrai0/m9-pca/notebook

library(dplyr)
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

# Get disorder labels for each sample
labels <- metadata$refinebio_subject

# Do the PCA
pc <- prcomp(gene_expression_t , scale=TRUE)

# Extract the first four principal components
pca_scores <- as.data.frame(pc$x[, 1:4])

# Name the columns for easy reference
colnames(pca_scores) <- c("PC1", "PC2", "PC3", "PC4")

# Combine PCA scores with disorder labels
pca_data <- cbind(pca_scores, disorder = labels)

# Extract the proportion of variance explained by the first two PCs
var_percentages <- summary(pc)$importance[2, 1:4] * 100  # in percentage

# Plot PC1 vs PC2
plot_pc1_pc2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = disorder)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Gene Expression Data (PC1 vs PC2)",
       x = paste("Principal Component 1 (", sprintf("%.2f", var_percentages[1]), "%)", sep = ""),
       y = paste("Principal Component 2 (", sprintf("%.2f", var_percentages[2]), "%)", sep = ""),
       color = "disorder") +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Plot PC1 vs PC3
plot_pc1_pc3 <- ggplot(pca_data, aes(x = PC1, y = PC3, color = disorder)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Gene Expression Data (PC1 vs PC3)",
       x = paste("Principal Component 1 (", sprintf("%.2f", var_percentages[1]), "%)", sep = ""),
       y = paste("Principal Component 3 (", sprintf("%.2f", var_percentages[3]), "%)", sep = ""),
       color = "disorder") +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Plot PC1 vs PC4
plot_pc1_pc4 <- ggplot(pca_data, aes(x = PC1, y = PC4, color = disorder)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Gene Expression Data (PC1 vs PC4)",
       x = paste("Principal Component 1 (", sprintf("%.2f", var_percentages[1]), "%)", sep = ""),
       y = paste("Principal Component 4 (", sprintf("%.2f", var_percentages[4]), "%)", sep = ""),
       color = "disorder") +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Plot PC2 vs PC3
plot_pc2_pc3 <- ggplot(pca_data, aes(x = PC2, y = PC3, color = disorder)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Gene Expression Data (PC2 vs PC3)",
       x = paste("Principal Component 2 (", sprintf("%.2f", var_percentages[2]), "%)", sep = ""),
       y = paste("Principal Component 3 (", sprintf("%.2f", var_percentages[3]), "%)", sep = ""),
       color = "disorder") +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Plot PC2 vs PC4
plot_pc2_pc4 <- ggplot(pca_data, aes(x = PC2, y = PC4, color = disorder)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Gene Expression Data (PC2 vs PC4)",
       x = paste("Principal Component 2 (", sprintf("%.2f", var_percentages[2]), "%)", sep = ""),
       y = paste("Principal Component 4 (", sprintf("%.2f", var_percentages[4]), "%)", sep = ""),
       color = "disorder") +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Show the plots
plot_pc1_pc2
plot_pc1_pc3
plot_pc1_pc4
plot_pc2_pc3
plot_pc2_pc4

# Save the plots
ggsave("./plots/PCA/PCA-PC1-vs-PC2.png", plot_pc1_pc2, width = 16, height = 9)
ggsave("./plots/PCA/PCA-PC1-vs-PC3.png", plot_pc1_pc3, width = 16, height = 9)
ggsave("./plots/PCA/PCA-PC1-vs-PC4.png", plot_pc1_pc4, width = 16, height = 9)
ggsave("./plots/PCA/PCA-PC2-vs-PC3.png", plot_pc2_pc3, width = 16, height = 9)
ggsave("./plots/PCA/PCA-PC2-vs-PC4.png", plot_pc2_pc4, width = 16, height = 9)

# Summary of the PCA
summary(pc)

# Plot each PC against the proportion of variance explained

# Extract the proportion of variance explained
var_explained <- summary(pc)$importance["Proportion of Variance", ]

# Convert the list to a named vector
my_vector <- unlist(var_explained)

png("./plots/PCA/PCA-barplot.png")

# Create a bar plot for top 20 PCs
pc_barplot <- barplot(my_vector[1:20], main="Bar Plot of Top 20 Principle Components", xlab="Principle Components", ylab="Percentage of Variance Explained", col="blue")

dev.off()
