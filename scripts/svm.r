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

# All predictions
all_predictions <- data.frame()

# Create total confusion matrix
total_confusion_matrix <- matrix(0, nrow=2, ncol=2)

# ----------------------------------------------------------------------

## SVM with 5-fold cross validation

for(i in 1:num_folds) {

  # Get indices for training and testing datasets
  test_indices <- folds[[i]]
  train_indices <- unlist(folds[-i])

  # Split into training and testing datasets
  train <- gene_expression_t[train_indices,]
  test <- gene_expression_t[test_indices,]

  # Train SVM

  # LINEAR
  svm_model <- svm(disorder ~ ., data=train, kernel="linear")

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

  # Create a confusion matrix
    confusion_matrix <- table(test$disorder, predictions)
    print(confusion_matrix)
    total_confusion_matrix <- total_confusion_matrix + confusion_matrix

  # Print out predictions
  predictions_df <- as.data.frame(predictions)
  all_predictions <- rbind(all_predictions, predictions_df)

}

# Print out total confusion matrix
cat("Total Confusion Matrix:\n")
print(total_confusion_matrix)

# Calculate average AUC of all trials
average_auc <- mean(auc_values)

# Calculate precision, recall, and F1 for predicting positive class (schizophrenia)
positive_precision <- total_confusion_matrix[2,2] / (total_confusion_matrix[2,2] + total_confusion_matrix[1,2])
positive_recall <- total_confusion_matrix[2,2] / (total_confusion_matrix[2,2] + total_confusion_matrix[2,1])
positive_f1 <- 2 * (positive_precision * positive_recall) / (positive_precision + positive_recall)

# Calculate precision, recall, and F1 for predicting negative class (control)
negative_precision <- total_confusion_matrix[1,1] / (total_confusion_matrix[1,1] + total_confusion_matrix[2,1])
negative_recall <- total_confusion_matrix[1,1] / (total_confusion_matrix[1,1] + total_confusion_matrix[1,2])
negative_f1 <- 2 * (negative_precision * negative_recall) / (negative_precision + negative_recall)

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

# Arrange predictions by sample ID and write predictions to file
all_predictions <- all_predictions %>% arrange(rownames(all_predictions))
all_predictions <- all_predictions %>%
  tibble::rownames_to_column("Gene")
readr::write_tsv(all_predictions, "./results/svm/linear-svm-5000-predictions.tsv")


# # ----------------------------------------------------------------------

# # EXTRACT GENE SIGNATURES

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
# # Convert rownames to column
# gene_weights_sorted <- gene_weights_sorted %>%
#   tibble::rownames_to_column("Gene")

# # Plot the sorted weights
# plot(gene_weights_sorted$weight, ylab="Values", xlab="Index", main="Plot of Values", pch=16, cex=0.5)

# # # Filter for top X genes
# top_x_genes <- 200
# gene_weights_sorted <- gene_weights_sorted[1:top_x_genes,]
# # Get gene names
# gene_signatures <- gene_weights_sorted$Gene

# # Write gene signatures to file
# readr::write_tsv(as.data.frame(gene_signatures), "./results/svm/linear-svm-5000-genes-200-signatures.tsv")


# # ----------------------------------------------------------------------

# # Compute the area under the ROC curve
# decision_values <- predict(svm_model, gene_expression_t, decision.values=TRUE)
# attr_values <- attr(decision_values, "decision.values")
# attr_values <- as.vector(attr_values)
# roc_obj <- roc(response=gene_expression_t$disorder, predictor=attr_values)
# auc(roc_obj)
# plot(roc_obj, main="ROC Curve")



