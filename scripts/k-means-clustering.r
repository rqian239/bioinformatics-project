# Perform k-means clustering on gene expression data
# From https://github.com/SpencerPao/Data_Science/blob/main/KMeans/kmeans.R

# install.packages("factoextra")

# Load libraries
library(factoextra)
library(dplyr)

# Read in gene expression data
gene_expression <- readr::read_tsv("./data/SRP073813/SRP073813-HUGO-cleaned-top-5000.tsv")

# Transpose the data
gene_expression_t <- as.data.frame(t(gene_expression))

# Make the first row the column names
colnames(gene_expression_t) <- gene_expression_t[1,]
gene_expression_t <- gene_expression_t[-1,]

# Load in metadata
metadata <- readr::read_tsv("./data/SRP073813/metadata_SRP073813_filtered.tsv") %>% tibble::column_to_rownames("refinebio_accession_code")
# Check rownames match
all(rownames(metadata) == rownames(gene_expression_t))

# Get labels of samples
sample_labels <- metadata$refinebio_subject
table(sample_labels)


# Scale the data (this is not needed anymore since the data has been log-scaled and normalized)
gene_expression_t <- mutate_all(gene_expression_t, function(x) as.numeric(as.character(x)))
gene_expression_t_scaled <- gene_expression_t
# gene_expression_t_scaled <- as.data.frame(scale(gene_expression_t))

# Distance
gene_expression_t_dist <- dist(gene_expression_t_scaled)

# Calculate how many clusters we should have, using the elbow method
fviz_nbclust(gene_expression_t_scaled, kmeans, method = "wss") +
    labs(subtitle = "Elbow method")

# Perform k-means clustering
k <- 5
km.out <- kmeans(gene_expression_t_scaled, centers = k, nstart = 55)    # Go through each of the 55 possible starting points
print(km.out)

# Table of cluster and group
table(sample_labels, km.out$cluster)

# Visualize the clusters
km.clusters <- km.out$cluster
gene_expression_t_plot <- gene_expression_t_scaled
rownames(gene_expression_t_plot) <- paste(sample_labels, 1:dim(gene_expression_t_plot)[1], sep = "_")
fviz_cluster(list(data=gene_expression_t_plot, cluster=km.clusters))

# Write the cluster assignments to a file
cluster_assignment_results <- data.frame(km.out$cluster)
colnames(cluster_assignment_results)[1] <- "k_cluster"

# Check if metadata and cluster assignments match
all(rownames(metadata) == rownames(cluster_assignment_results))

# Add psychiatric disorder labels to the cluster assignments
cluster_assignment_results <- cbind(cluster_assignment_results, disorder_group = metadata$refinebio_subject)
cluster_assignment_results <- rownames_to_column(cluster_assignment_results, "refinebio_accession_code")

# Write the cluster assignments to a file
readr::write_tsv(cluster_assignment_results, "./results/clustering/k-means-cluster-results.tsv")