# Load in clustering results
cluster_df <- readr::read_tsv("./results/clustering/k-means-cluster-results.tsv")

table(cluster_df$disorder_group, cluster_df$k_cluster)

# Perform chi-squared test
test_result <- chisq.test(cluster_df$disorder_group, cluster_df$k_cluster)

# adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
