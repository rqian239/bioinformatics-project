# install.packages("ggalluvial")

library(ggalluvial)
library(dplyr)

# Load in cluster data
k_means_clusters <- readr::read_tsv("./results/clustering/k-means-cluster-results.tsv")

# Combine data into one dataframe to plot
clustering_results <- k_means_clusters # Combine clustering results here


# Create alluvial plot
ggplot(data = clustering_results,
    aes(axis1 = disorder_group, axis2 = k_cluster)) +
    scale_x_discrete(limits = c("disorder_group", "k_cluster"), expand = c(.2, .05)) + 
    xlab("Grouping") +
    geom_alluvium(aes(fill = disorder_group)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) + 
    theme_minimal() + 
    ggtitle("Clustering Memberships of Different Methods")