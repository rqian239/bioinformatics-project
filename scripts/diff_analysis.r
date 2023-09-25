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

