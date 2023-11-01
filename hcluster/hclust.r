library(data.table)
library(ggplot2)

gene_expression <- read.csv(file = "./data/SRP073813/gene_expression_schizophrenia.csv", row.names = 1)
gene_metadata <- read.csv(file = "./data//SRP073813/gene_metadata.csv")

gene_metadata <- gene_metadata[-1, ]

gene_expression <- as.data.frame(gene_expression)
gene_expression <- gene_expression[order(apply(gene_expression, 1, var), decreasing = TRUE), ]

for (i in 1:5)
{
    x <- 10000

    if (i == 2) {
        x <- 1000
    } else if (i == 3) {
        x <- 100
    } else if (i == 4) {
        x <- 10
    } else if (i ==5) {
        x <- 5000
    }

    disimilarity <- dist(gene_expression[1:x], method = "euclidean")
    hcl <- hclust(disimilarity, method = "ward.D2")

    filename <- sprintf("./hcluster/Hierarchical_Cluster_%d.pdf", x)
    pdf(filename, width = 16, height = 9)
    plot(hcl, main = filename)

    dev.off()
}

# k_clusters <- as.data.frame(cutree(hcl, k = 5))
# k_clusters$subject_id <- row.names(k_clusters)
# condition <- data.frame()

# for (subj_id in k_clusters$subject_id)
# {
#     disorder_group <- as.character((gene_metadata[grep(subj_id, gene_metadata$refinebio_accession_code), ])$refinebio_disease)
#     condition <- rbind(condition, disorder_group)
# }
# row.names(condition) <- k_clusters$subject_id
# colnames(condition)[1] <- "disorder_group"

# k_clusters <- cbind(k_clusters, condition)

# colnames(k_clusters)[1] <- "hierarchical_cluster"
# k_clusters <- k_clusters[-2]
# write.csv(k_clusters, file = "./hcluster/h_cluster.csv")
