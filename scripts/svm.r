# Support Vector Machine

# install.packages("e1071")
# install.packages("pROC")

# Load libraries
library(e1071)
library(dplyr)
library(pROC)

data_filepath <- "./data/SRP073813/SRP073813-HUGO-cleaned-top-5000.tsv"

# Read in gene expression data
gene_expression <- readr::read_tsv(data_filepath)

# Transpose the data
gene_expression_t <- as.data.frame(t(gene_expression))

# Make the first row the column names
colnames(gene_expression_t) <- gene_expression_t[1,]
gene_expression_t <- gene_expression_t[-1,]

# Convert all numbers into a numeric
gene_expression_t <- mutate_all(gene_expression_t, function(x) as.numeric(as.character(x)))

# Load in metadata
metadata <- readr::read_tsv("./data/SRP073813/metadata_SRP073813_filtered.tsv") %>% tibble::column_to_rownames("refinebio_accession_code")
# Check rownames match
all(rownames(metadata) == rownames(gene_expression_t))

# Add column for psychiatric disorder into gene expression data
gene_expression_t$disorder <- as.factor(metadata$refinebio_subject)


## SVM

# Linear SVM model with 5-fold cross validation
svm_model <- svm(disorder ~ ., data=gene_expression_t, cross=5, kernel="linear")


# Extract gene signatures
weights <- t(svm_model$coefs) %*% svm_model$SV
# Associate weights with gene names
gene_weights <- as.data.frame(t(setNames(weights, rownames(gene_expression_t))))
# Rename column in gene_weights
colnames(gene_weights)[colnames(gene_weights) == "V1"] <- "weight"
# Absolute value all weights
gene_weights$weight <- abs(gene_weights$weight)
# Sort by weight
gene_weights_sorted <- gene_weights %>% arrange(desc(weight))


# Compute the area under the ROC curve
decision_values <- predict(svm_model, gene_expression_t, decision.values=TRUE)
attr_values <- attr(decision_values, "decision.values")
attr_values <- as.vector(attr_values)
roc_obj <- roc(response=gene_expression_t$disorder, predictor=attr_values)
auc(roc_obj)
plot(roc_obj, main="ROC Curve")



