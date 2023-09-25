# Script from https://alexslemonade.github.io/refinebio-examples/03-rnaseq/differential-expression_rnaseq_01.html

# Install packages if not already installed
if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}

library(DESeq2)
library(ggplot2)
library(magrittr)
library(dplyr)

# Set seed for DESeq2::plotCounts() function
set.seed(12345)

# Read in metadata TSV file
metadata_file <- "./data/SRP073813/metadata_SRP073813.tsv"
metadata <- readr::read_tsv(metadata_file)

# Read in gene expression data
data_file <- "./data/SRP073813/SRP073813-HUGO-cleaned.tsv"
expression_df <- readr::read_tsv(data_file) %>% 
  mutate(first_mapped_hugo = make.unique(as.character(first_mapped_hugo))) %>%
  tibble::column_to_rownames("first_mapped_hugo")

# Note: Duplicate gene names were found in the data set. Duplicates were renamed with make.unique(), a function that appends a suffix to each duplicate.

# DUPLICATE GENE NAMES: 'ARL17A', 'ASPRV1', 'BAZ2B', 'BOLA2', 'C2-AS1', 'C4orf36', 'CHN2', 'CLEC4A', 
# 'CLN3', 'COMMD7', 'CRHR1', 'CT47A11', 'CYB561D2', 'DDT', 'DEFB105A', 'DEFB4A', 'DET1', 'DGCR6', 
# 'DNAJC9-AS1', 'ELFN2', 'FAM163A', 'FAM50B', 'GGT1', 'GOLGA6L4', 'GOLGA8M', 'HLA-DQA1', 'HNRNPA1L2', 'HOXC5', 
# 'IFNAR2', 'ITFG2-AS1', 'KBTBD11-OT1', 'KYAT1', 'LGALS7', 'LINC00486', 'LINC01115', 'LINC01238', 'LINC02203', 'LINC03021', 
# 'LINC03025', 'LOC105377488', 'LRFN2', 'LRRC57', 'MAFIP', 'MATR3', 'MIEF1', 'MKKS', 'MPV17L', 'MRPL23', 'NDST2', 'OR1F1', 
# 'PDE4C', 'PDE8B', 'PINX1', 'POLA2', 'PRH1', 'PRODH', 'PRORP', 'PTGR2', 'PTP4A1', 'RAET1E', 'RAET1E-AS1', 'RBM27', 'RGS5', 
# 'RN7SK', 'SCHIP1', 'SEC16B', 'SERF1A', 'SFTA3', 'SLC35D2', 'SMN1', 'SPATA13', 'SPICE1', 'SRSF1', 'ST6GALNAC6', 'STXBP2', 
# 'TAB2', 'TBC1D3', 'TBCE', 'TMSB15B', 'TUT1', 'USP17L8', 'ZNF177', 'ZNF280D', 'ZNF544', 'ZNF724', 'ZNF844'

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# Print the values we need for comparison
head(metadata$refinebio_title)

# Label the psychiatric disorder using the metadata

# Comparing all 4 groups here?
# metadata <- metadata %>%
#   dplyr::mutate(psy_disorder = dplyr::case_when(
#     stringr::str_detect(refinebio_title, "_C_") ~ "Control",
#     stringr::str_detect(refinebio_title, "_M_") ~ "Major Depression",
#     stringr::str_detect(refinebio_title, "_B_") ~ "Bipolar Disorder",
#     stringr::str_detect(refinebio_title, "_S_") ~ "Schizophrenia"
#   ))

# Comparing control vs any disorder
metadata <- metadata %>%
  dplyr::mutate(psy_disorder = dplyr::case_when(
    stringr::str_detect(refinebio_title, "_C_") ~ "Control",
    TRUE ~ "Disorder"
  ))

# Label brain region using the metadata
metadata <- metadata %>%
  dplyr::mutate(brain_region = dplyr::case_when(
    stringr::str_detect(refinebio_title, "AnCg") ~ "Anterior Cingulate Cortex",
    stringr::str_detect(refinebio_title, "nAcc") ~ "Nucleus Accumbens",
    stringr::str_detect(refinebio_title, "DLPFC") ~ "Dorsolateral Prefrontal Cortex"
  ))

# Let's take a look at the original metadata column's info
# and our new columns
dplyr::select(metadata, refinebio_title, psy_disorder, brain_region)

# Print out a preview of `psy_disorder`
str(metadata$psy_disorder)


# Make mutation_status a factor and set the levels appropriately

# Comparing all 4 groups here?
# metadata <- metadata %>%
#   dplyr::mutate(
#     # Here we define the values our factor variable can have and their order.
#     psy_disorder = factor(psy_disorder, levels = c("Control", "Major Depression", "Bipolar Disorder", "Schizophrenia"))
#   )

# Comparing control vs any disorder
metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    psy_disorder = factor(psy_disorder, levels = c("Control", "Disorder"))
  )

# Check the levels
levels(metadata$psy_disorder)

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff

min_counts <- 10

# Filter the data
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= min_counts)

# Creating a DESeq2DataSet

# round all expression counts
gene_matrix <- round(filtered_expression_df)

# 0 COUNTS FOR EACH ROW, REPLACE WITH 1
# I get an error running DESeq2: Error: Every gene contains at least one zero, cannot compute log geometric means
# A quick solution is to add 1 for every zero count in the data set. This will be 0 once we take the log.
gene_matrix[gene_matrix == 0] <- 1

# Create a DESeq2DataSet object
ddset <- DESeqDataSetFromMatrix(
  # Here we supply non-normalized count data
  countData = gene_matrix,
  # Supply the `colData` with our metadata data frame
  colData = metadata,
  # Supply our experimental variable to `design`
  design = ~psy_disorder
)

# Run the differential expression analysis
deseq_object <- DESeq(ddset)

# Extract the results
deseq_results <- results(deseq_object)

# Obtain shrunken log fold change estimates based on negative binomial distribution
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq()   NOTE: this is the coefficient is how many groups you have
  res = deseq_results # The original DESeq2 results table
)

# Print some of the results
head(deseq_results)

# Save results into dataframe
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

head(deseq_df)


# Plot one gene to show difference
plotCounts(ddset, gene = "SERPINA3", intgroup = "psy_disorder")


# Write results to tsv file
readr::write_tsv(
deseq_df,
"./results/SRP073813_diff_expr_results.tsv"
)


# CREATE VOLCANO PLOT

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

# Print out plot here
volcano_plot


# Extract differentially expressed genes


  