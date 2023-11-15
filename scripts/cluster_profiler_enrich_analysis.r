# From https://github.com/mousepixels/sanbomics_scripts/blob/main/GO_in_R.Rmd

# Run Gene Ontology (GO) enrichment analysis using clusterProfiler

# BiocManager::install("clusterProfiler")
# BiocManager::install("AnnotationDbi")

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Get the significantly differentially expressed genes
diff_expressed_genes <- readr::read_tsv("./results/diff_expression/SRP073813_diff_expressed_genes.tsv")

# Get gene names
genes <- diff_expressed_genes$Gene

print(genes)

# Run enrichGO
GO_results <- enrichGO(gene = genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")

# Convert results to data.frame
GO_results_df <- as.data.frame(GO_results)

# Order results by p-value
GO_results_df <- GO_results_df %>% arrange(pvalue)

write.table(GO_results_df, file = "./results/diff_expression/clust_profiler_GO_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# Plot the results
fit <- plot(barplot(GO_results, showCategory = 15))

png("./plots/diff_expression/cluster_profiler.png", res = 250, width = 1400, height = 1800)
print(fit)
dev.off()

