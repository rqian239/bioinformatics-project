# Support Vector Machine

# install.packages("e1071")
# install.packages("pROC")
# install.packages("caret")

# Load libraries
library(e1071)
library(dplyr)
library(pROC)
library(caret)

# DATA FORMATTING

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

# ----------------------------------------------------------------------


# Set up for 5-fold cross validation
num_folds <- 5

# Create folds
folds <- createFolds(gene_expression_t$disorder, k = num_folds)

# ROC list for each fold
roc_list <- list()

# AUC values for each fold
auc_values <- numeric(num_folds)

# ----------------------------------------------------------------------

## SVM

for(i in 1:num_folds) {

  # Get indices for training and testing datasets
  test_indices <- folds[[i]]
  train_indices <- unlist(folds[-i])

  # Split into training and testing datasets
  train <- gene_expression_t[train_indices,]
  test <- gene_expression_t[test_indices,]

  # Train SVM

#   # LINEAR
#   svm_model <- svm(disorder ~ ., data=train, kernel="linear")

#   # POLYNOMIAL
#     svm_model <- svm(disorder ~ ., data=train, kernel="polynomial")

#    # RADIAL (RBF)
#     svm_model <- svm(disorder ~ ., data=train, kernel="radial")

#   # SIGMOID
#     svm_model <- svm(disorder ~ ., data=train, kernel="sigmoid")

  # Predict on test set
  predictions <- predict(svm_model, test, decision.values=TRUE)

  # Extract decision values (scores) for the ROC curve
  decision_values <- attr(predictions, "decision.values")[,1]

  # Create ROC curve
  roc_obj <- roc(test$disorder, decision_values)

  # Store the ROC object in the list
  roc_list[[i]] <- roc_obj

  # Store the AUC value
  auc_values[i] <- auc(roc_obj)

#   # Plot ROC curve
#   plot(roc_obj, main=sprintf("ROC Curve for Fold %d", i), col=i)

  # Print AUC
  cat(sprintf("AUC for Fold %d: %f\n", i, auc(roc_obj)))

}

average_auc <- mean(auc_values)

# Plot all ROC curves
if (length(folds) > 1) {
  plot(roc_list[[1]], main="ROC Curves for All Folds")
  for (i in 2:length(folds)) {
    lines(roc_list[[i]], col=i)
  }
}

# Add a legend
legend_labels <- sprintf("Fold %d (AUC = %.2f)", 1:length(folds), auc_values)
legend("bottomright", legend=legend_labels, col=1:length(folds), lty=1)




# # Linear SVM model with 5-fold cross validation
# svm_model <- svm(disorder ~ ., data=gene_expression_t, cross=5, kernel="linear")


# # Extract gene signatures
# weights <- t(svm_model$coefs) %*% svm_model$SV
# # Associate weights with gene names
# gene_weights <- as.data.frame(t(setNames(weights, rownames(gene_expression_t))))
# # Rename column in gene_weights
# colnames(gene_weights)[colnames(gene_weights) == "V1"] <- "weight"
# # Absolute value all weights
# gene_weights$weight <- abs(gene_weights$weight)
# # Sort by weight
# gene_weights_sorted <- gene_weights %>% arrange(desc(weight))


# # Compute the area under the ROC curve
# decision_values <- predict(svm_model, gene_expression_t, decision.values=TRUE)
# attr_values <- attr(decision_values, "decision.values")
# attr_values <- as.vector(attr_values)
# roc_obj <- roc(response=gene_expression_t$disorder, predictor=attr_values)
# auc(roc_obj)
# plot(roc_obj, main="ROC Curve")



