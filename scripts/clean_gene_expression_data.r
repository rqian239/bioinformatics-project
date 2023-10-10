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

# Filter metadata to only include samples from the Anterior Cingulate Cortex (AnCg) from Control and Schizophrenia patients
selected_samples <- metadata %>%
  dplyr::filter(refinebio_subject %in% c("ancg_control", "ancg_schizophrenia")) %>%
  dplyr::select(refinebio_accession_code)

# Filter gene expression data to only include samples from the Anterior Cingulate Cortex (AnCg) from Control and Schizophrenia patients
final_gene_expression_data <- gene_expression_data_cleaned %>%
    dplyr::select(first_mapped_hugo, all_of(selected_samples$refinebio_accession_code))

# Write the cleaned/filtered data to a new file
write.table(final_gene_expression_data, "SRP073813-HUGO-cleaned.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
