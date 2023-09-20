library(magrittr)
library(org.Hs.eg.db)
metadata <- readr::read_tsv("./SRP073813/metadata_SRP073813.tsv")

expression_df <- readr::read_tsv("./SRP073813/SRP073813.tsv") %>%
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

mapped_list <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "list"
)

mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "HUGO") %>%
  # enframe() makes a `list` column; we will simplify it with unnest()
  # This will result in one row of our data frame per list item
  tidyr::unnest(cols = HUGO)


collapsed_mapped_df <- mapped_df %>%
  # Group by Ensembl IDs
  dplyr::group_by(Ensembl) %>%
  # Collapse the HUGO SYMBOLS `mapped_df` into one column named `all_hugo_ids`
  dplyr::summarize(all_hugo_ids = paste(HUGO, collapse = ";"))

collapsed_mapped_df %>%
  # Filter `collapsed_mapped_df` to include only the rows where
  # `all_hugo_ids` values include the ";" character --
  # these are the rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_hugo_ids, ";")) %>%
  # We only need a preview here
  head()

final_mapped_df <- data.frame(
  "first_mapped_entrez_id" = mapIds(
    org.Hs.eg.db, # Replace with annotation package for your organism
    keys = expression_df$Gene,
    keytype = "ENSEMBL", # Replace with the gene identifiers used in your data
    column = "SYMBOL", # The type of gene identifiers you would like to map to
    multiVals = "first" # Keep only the first mapped value for each Ensembl ID
  )
) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Add the multiple mappings data from `collapsed_mapped_df` using Ensembl IDs
  dplyr::inner_join(collapsed_mapped_df, by = "Ensembl") %>%
  # Now let's add on the rest of the expression data
  dplyr::inner_join(expression_df, by = c("Ensembl" = "Gene"))


readr::write_tsv(final_mapped_df, file.path(
  "./",
  "SRP073813-HUGO.tsv" # Replace with a relevant output file name
))

