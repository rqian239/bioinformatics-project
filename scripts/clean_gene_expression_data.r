# install.packages("data.table")
library(data.table)

# Load in the gene expression data that has already been mapped to HUGO
gene_expression_data <- fread("./data/SRP073813/SRP073813-HUGO.tsv", header = TRUE, sep = "\t")

# Remove the `Ensembl` and `all_hugo_ids` columns
columns_to_remove <- c("Ensembl", "all_hugo_ids")
gene_expression_data[, (columns_to_remove) := NULL]

# Remove rows with "NA" HUGO IDs
gene_expression_data_cleaned <- gene_expression_data[!is.na(gene_expression_data$first_mapped_hugo)]

# Write the cleaned data to a new file
write.table(gene_expression_data_cleaned, "SRP073813-HUGO-cleaned.tsv", sep = "\t", quote = FALSE, row.names = FALSE)