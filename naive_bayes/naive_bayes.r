library(data.table)
library(ggplot2)
library(dplyr)
library(parsnip)
library(rsample)
library(discrim)
library(recipes)
library(tune)
library(naivebayes)
library(tidymodels)
library(vip)


gene_no <- 10000
gene_metadata <- read.csv(file = "./data//SRP073813/gene_metadata.csv")
gene_metadata <- gene_metadata[-1, ]
gene_metadata <- gene_metadata %>% filter(refinebio_subject %in% c("ancg_control", "ancg_schizophrenia"))
gene_metadata <- gene_metadata %>% select(c("refinebio_accession_code", "refinebio_disease"))


gene_expression <- read.csv(file = "./data/SRP073813/gene_expression_schizophrenia.csv", row.names = 1)

# Select the columns from the data frame in the order specified by the new vector
gene_expression <- gene_expression[order(apply(gene_expression, 1, var), decreasing = TRUE), ]

gene_expression <- gene_expression[1:gene_no]
gene_expression$refinebio_accession_code <- row.names(gene_expression)

gene_expression <- gene_expression %>% merge(gene_metadata, by = "refinebio_accession_code")
colnames(gene_expression)[gene_no+2] <- "class"

rownames(gene_expression) <- gene_expression$refinebio_accession_code
# Remove the first column of the gene_expression data frame
gene_expression <- gene_expression[-1]


set.seed(123)

gene_recipe <- recipe(class ~ ., data = gene_expression)

cross_validation <- vfold_cv(gene_expression, v = 5)

# Creating the model

nb_model <- naive_Bayes() %>% 
    set_mode("classification") %>%
    set_engine("naivebayes")

extract <- function(x) {
    extract_fit_parsnip(x)
}



results <- fit_resamples(nb_model, gene_recipe, cross_validation, control = control_resamples(save_pred = TRUE, extract = extract))

metrics <- collect_metrics(results, summarize = FALSE)


extract_results <- results$.extracts[[1]]$.extracts[[1]]
