library(data.table)
library(ggplot2)
library(tidymodels)
library(discrim)

gene_expression <- read.csv(file = "./data/SRP073813/gene_expression_schizophrenia.csv", row.names = 1)

gene_expression$refinebio_accession_code <- row.names(gene_expression)


gene_metadata <- read.csv(file = "./data//SRP073813/gene_metadata.csv")
gene_metadata <- gene_metadata[-1, ]
gene_metadata <- gene_metadata %>% filter(refinebio_subject %in% c("ancg_control", "ancg_schizophrenia"))
gene_metadata <- gene_metadata %>% select(c("refinebio_accession_code", "refinebio_disease"))

gene_expression <- gene_expression %>% merge(gene_metadata, by = "refinebio_accession_code")


variances <- var(gene_expression)

# Sort the columns in the data frame by their variance in descending order
gene_expression <- gene_expression %>% select(order(variances, decreasing = TRUE))



# Split your data into training and testing sets
# initial_split <- initial_split(gene_expression, prop = 3/4)
# training_data <- training(initial_split)
# testing_data <- testing(initial_split)