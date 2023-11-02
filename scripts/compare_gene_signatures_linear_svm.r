# Read all gene signatures

gene_sig_10000 <- readr::read_tsv("./results/svm/linear-svm-10000-gene-signatures.tsv")
gene_sig_5000 <- readr::read_tsv("./results/svm/linear-svm-5000-genes-200-signatures.tsv")
gene_sig_1000 <- readr::read_tsv("./results/svm/linear-svm-1000-gene-signatures.tsv")
gene_sig_100 <- readr::read_tsv("./results/svm/linear-svm-100-gene-signatures.tsv")
gene_sig_10 <- readr::read_tsv("./results/svm/linear-svm-10-gene-signatures.tsv")

# Compare gene signatures with 5000 trial

common_gene_sigs_5000_10000 <- intersect(gene_sig_5000$gene_signatures, gene_sig_10000$gene_signatures)
print(length(common_gene_sigs_5000_10000))

common_gene_sigs_5000_1000 <- intersect(gene_sig_5000$gene_signatures, gene_sig_1000$gene_signatures)
print(length(common_gene_sigs_5000_1000))

common_gene_sigs_5000_100 <- intersect(gene_sig_5000$gene_signatures, gene_sig_100$gene_signatures)
print(length(common_gene_sigs_5000_100))

common_gene_sigs_5000_10 <- intersect(gene_sig_5000$gene_signatures, gene_sig_10$gene_signatures)
print(length(common_gene_sigs_5000_10))

# Compare gene signatures with 1000 trial

common_gene_sigs_1000_10000 <- intersect(gene_sig_1000$gene_signatures, gene_sig_10000$gene_signatures)
print(length(common_gene_sigs_1000_10000))

common_gene_sigs_1000_100 <- intersect(gene_sig_1000$gene_signatures, gene_sig_100$gene_signatures)
print(length(common_gene_sigs_1000_100))

common_gene_sigs_1000_10 <- intersect(gene_sig_1000$gene_signatures, gene_sig_10$gene_signatures)
print(length(common_gene_sigs_1000_10))

# Compare gene signatures with 100 trial

common_gene_sigs_100_10000 <- intersect(gene_sig_100$gene_signatures, gene_sig_10000$gene_signatures)
print(length(common_gene_sigs_100_10000))

common_gene_sigs_100_10 <- intersect(gene_sig_100$gene_signatures, gene_sig_10$gene_signatures)
print(length(common_gene_sigs_100_10))

# Compare gene signatures with 10 trial

common_gene_sigs_10_10000 <- intersect(gene_sig_10$gene_signatures, gene_sig_10000$gene_signatures)
print(length(common_gene_sigs_10_10000))

# Compare gene signatures with 10000 trial

