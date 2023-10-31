# install.packages("ggalluvial")

library(ggalluvial)
library(dplyr)

# Load in cluster data
clustering_results <- readr::read_tsv("./results/clustering/all-cluster-results.tsv")

clustering_results$k_cluster <- paste("k-", clustering_results$k_cluster, sep = "")
clustering_results$hierarchical_cluster <- paste("h-", clustering_results$hierarchical_cluster, sep = "")
clustering_results$pam_cluster <- paste("p-", clustering_results$pam_cluster, sep = "")

# Create alluvial plot
ggplot(data = clustering_results,
    aes(axis1 = disorder_group, axis2 = k_cluster, axis3 = hierarchical_cluster, axis4 = pam_cluster)) +
    scale_x_discrete(limits = c("disorder_group", "k_cluster", "hierarchical_cluster", "pam_cluster"), expand = c(.2, .05)) + 
    xlab("Grouping") +
    geom_alluvium(aes(fill = disorder_group)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) + 
    theme_minimal() + 
    ggtitle("Clustering Memberships of Different Methods")