# From https://www.biostars.org/p/350710/

# Use topGO to perform GO enrichment analysis

library("topGO")
library("org.Hs.eg.db")

# Get the significantly differentially expressed genes
diff_expressed_genes <- readr::read_tsv("./results/SRP073813_diff_expressed_genes.tsv")

# Get gene names
genes <- diff_expressed_genes$Gene


