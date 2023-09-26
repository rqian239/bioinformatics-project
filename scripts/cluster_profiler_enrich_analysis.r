# From https://github.com/mousepixels/sanbomics_scripts/blob/main/GO_in_R.Rmd

# Run Gene Ontology (GO) enrichment analysis using clusterProfiler

# BiocManager::install("clusterProfiler")
# BiocManager::install("AnnotationDbi")

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Get the significantly differentially expressed genes
diff_expressed_genes <- readr::read_tsv("./results/SRP073813_diff_expressed_genes.tsv")

# Get gene names
genes <- diff_expressed_genes$Gene

print(genes)

# Run enrichGO
GO_results <- enrichGO(gene = genes, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")

# Convert results to data.frame
as.data.frame(GO_results)

# Plot the results
fit <- plot(barplot(GO_results, showCategory = 15))

png("./plots/cluster_profiler.png", res = 250, width = 1400, height = 1800)
print(fit)
dev.off()

