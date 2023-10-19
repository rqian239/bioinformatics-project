# This script will combine all of the cluster data into one file

library(dplyr)

h_clust_results <- readr::read_tsv("./results/clustering/h-cluster-results.tsv")
k_clust_results <- readr::read_tsv("./results/clustering/k-means-cluster-results.tsv")
pam_clust_results <- readr::read_tsv("results/clustering/pam-results.tsv")

# # Check if sample order is the same
# have_identical_sample_IDs <- identical(h_clust_results$refinebio_accession_code, k_clust_results$refinebio_accession_code) && identical(h_clust_results$refinebio_accession_code, pam_clust_results$refinebio_accession_code)
# print(have_identical_sample_IDs)

# Rename pam cluster column
pam_clust_results <- pam_clust_results %>% rename(pam_cluster = pamx5_t.clustering)

# Merge all cluster results into one dataframe
all_cluster_results <- k_clust_results %>% 
    inner_join(h_clust_results, by = "refinebio_accession_code") %>% 
    inner_join(pam_clust_results, by = "refinebio_accession_code") %>% 
    select(refinebio_accession_code, disorder_group, k_cluster, hierarchical_cluster, pam_cluster)

# Write the cluster assignments to a file
readr::write_tsv(all_cluster_results, "./results/clustering/all-cluster-results.tsv")