# Load in clustering results
cluster_df <- read.csv("./hcluster/h_cluster.csv")

# Perform chi-squared test
chi <- chisq.test(cluster_df$disorder_group, cluster_df$hierarchical_cluster)

# adjusted_p_values <- p.adjust(chi$p.value, method = "bonferroni")

table(cluster_df$disorder_group, cluster_df$hierarchical_cluster)
