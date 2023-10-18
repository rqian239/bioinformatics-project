library(data.table)
library(ggplot2)

gene_expression <- read.csv(file = "./data/SRP073813/gene_expression_schizophrenia.csv")
gene_metadata <- read.csv(file = "./data/SRP073813/metadata_schizophrenia.tsv", sep = "\t")

# gene_metadata <- gene_metadata[-1,]
# gene_expression <- as.data.frame(t(gene_expression))

# names(gene_expression) <- gene_expression[1,]
# gene_expression <- gene_expression[-1,]
# gene_expression <- sapply(gene_expression, as.numeric)
# gene_expression <- as.data.frame(gene_expression)


# gene_expression$condition <- gene_metadata$refinebio_disease

gene_expression <- gene_expression[order(apply(gene_expression, 1, var), decreasing = TRUE), ]

for (i in 1:4)
{
    x <- 10000

    if (i == 2) {
        x <- 1000
    } else if (i == 3) {
        x <- 100
    } else if (i == 4) {
        x <- 10
    }

    disimilarity <- dist(gene_expression[1:x], method = "euclidean")
    hcl <- hclust(disimilarity, method = "ward.D2")

    filename <- sprintf("./cluster/Hierarchical_Cluster_%d.pdf", x)
    pdf(filename, width = 16, height = 9)
    plot(hcl, main = filename)

    dev.off()
}
