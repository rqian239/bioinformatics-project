library(gprofiler2)

genes <- read.csv(file = "./scripts/diff_expressed_genes.tsv", sep = "\t")

gostres <- gost(query = genes$Gene, organism = "hsapiens")

plot <- gostplot(gostres, interactive = FALSE)

publish_gostplot(plot, filename = "g_Profiler_plot.pdf")

publish_gostplot(plot, highlight_terms = gostres$result, width = 15, height = 30, filename = "g_Profiler_table.pdf")
