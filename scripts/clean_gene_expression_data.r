# This script will remove the 'EnsEMBL' and 'all_hugo_ids' columns from the gene expression data, and remove rows with "NA" HUGO IDs.
# This script will also remove any samples that were not collected from the Anterior Cingulate Cortex (AnCg) from Control and Schizophrenia patients.
# For our project, we are verifying the results of the original study to compare differences between Control and Schizophrenia patients in the AnCg.

# install.packages("data.table")
library(data.table)

# Load in the gene expression data that has already been mapped to HUGO
gene_expression_data <- fread("./data/SRP073813/SRP073813-HUGO.tsv", header = TRUE, sep = "\t")

# Remove the `Ensembl` and `all_hugo_ids` columns
columns_to_remove <- c("Ensembl", "all_hugo_ids")
gene_expression_data[, (columns_to_remove) := NULL]

# Remove rows with "NA" HUGO IDs
gene_expression_data_cleaned <- gene_expression_data[!is.na(gene_expression_data$first_mapped_hugo)]

# Now we will filter the data so that only samples collected from the Anterior Cingulate Cortex (ACC) from Control and Schizophrenia patients are included.

# Read in metadata TSV file
metadata_file <- "./data/SRP073813/metadata_SRP073813.tsv"
metadata <- readr::read_tsv(metadata_file)

# Label psychiatric disorder
metadata <- metadata %>%
  dplyr::mutate(psy_disorder = dplyr::case_when(
    stringr::str_detect(refinebio_title, "_C_") ~ "Control",
    stringr::str_detect(refinebio_title, "_M_") ~ "Major Depression",
    stringr::str_detect(refinebio_title, "_B_") ~ "Bipolar Disorder",
    stringr::str_detect(refinebio_title, "_S_") ~ "Schizophrenia"
  ))

# Label brain region
metadata <- metadata %>%
  dplyr::mutate(brain_region = dplyr::case_when(
    stringr::str_detect(refinebio_title, "AnCg") ~ "Anterior Cingulate Cortex",
    stringr::str_detect(refinebio_title, "nAcc") ~ "Nucleus Accumbens",
    stringr::str_detect(refinebio_title, "DLPFC") ~ "Dorsolateral Prefrontal Cortex"
  ))


# Write the cleaned data to a new file
write.table(gene_expression_data_cleaned, "SRP073813-HUGO-cleaned.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
