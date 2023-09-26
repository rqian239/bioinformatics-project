gene_exp <- read.csv("./data/SRP073813/SRP073813_HUGO.csv")
gene_exp <- gene_exp[(gene_exp$hugo_id != ""),]

gene_exp[,'X']
write.csv(gene_exp, file = "gene_expression.csv")