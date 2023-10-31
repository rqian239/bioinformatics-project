library(cluster)
library(dplyr)
library(ComplexHeatmap)
library(factoextra)

# generating pam clusters

df <- readr::read_tsv("./SRP073813-HUGO-cleaned-top-5000.tsv")

df_transpose <- as.data.frame(t(df))
colnames(df_transpose) <- df_transpose[1,]
df_transpose <- df_transpose[-1,]

metadata <- readr::read_tsv("./metadata_SRP073813_filtered.tsv") %>% tibble::column_to_rownames("refinebio_accession_code")
all(rownames(metadata) == rownames(df_transpose))

df_t2 <- df_transpose
rownames(df_t2) <- paste(as.list(metadata$refinebio_subject), 1:dim(metadata)[1], sep = "_")
df_t2 <- mutate_all(df_t2, function(x) as.numeric(as.character(x)))

pamx5_t <- pam(df_transpose, 5)
pdf("./pam-cluster-plots/pam-5.pdf")
fviz_cluster(list(data = df_t2, cluster = pamx5_t$clustering))
dev.off()

pamx5_t_1000 <- pam(df_transpose[, 1:1000], 5)
pdf("./pam-cluster-plots/pam-5_1000.pdf")
fviz_cluster(list(data = df_t2, cluster = pamx5_t_1000$clustering))
dev.off()

pamx5_t_100 <- pam(df_transpose[, 1:100], 5)
pdf("./pam-cluster-plots/pam-5_100.pdf")
fviz_cluster(list(data = df_t2, cluster = pamx5_t_100$clustering))
dev.off()

pamx5_t_10 <- pam(df_transpose[, 1:10], 5)
pdf("./pam-cluster-plots/pam-5_10.pdf")
fviz_cluster(list(data = df_t2, cluster = pamx5_t_10$clustering))
dev.off()

pamx7_t <- pam(df_transpose, 7)
pdf("./pam-cluster-plots/pam-7.pdf")
fviz_cluster(list(data = df_t2, cluster = pamx7_t$clustering))
dev.off()

pamx3_t <- pam(df_transpose, 3)
pdf("./pam-cluster-plots/pam-3.pdf")
fviz_cluster(list(data = df_t2, cluster = pamx3_t$clustering))
dev.off()

res <- data.frame(pamx5_t$clustering)
res <- cbind(res, disorder_group = metadata$refinebio_subject)
res <- tibble::rownames_to_column(res, "refinebio_accession_code")


readr::write_tsv(res, "./pam-results.tsv")


# heatmap
df <- readr::read_tsv("./SRP073813-HUGO-cleaned-top-5000.tsv") %>% tibble::column_to_rownames("first_mapped_hugo")

df_z <- t(apply(df, 1, scale))
colnames(df_z) <- colnames(df)

clustering_metadata <- readr::read_tsv("./pam-results.tsv") %>% tibble::column_to_rownames("refinebio_accession_code")
clustering_metadata <- clustering_metadata %>% arrange(pamx5_t.clustering)

grouping_order <- rownames(clustering_metadata)
heatmap_df <- df_z[, grouping_order]

anno <- HeatmapAnnotation(df = clustering_metadata, which = "col", 
    col = list(pamx5_t.clustering = c("1" = "blue", "2" = "red", "3" = "yellow", "4" = "green", "5" = "purple"), 
    disorder_group = c("ancg_control" = "lightcoral", "ancg_schizophrenia" = "cyan")))

pdf("heatmap.pdf", height=10, width=10)
Heatmap(heatmap_df, cluster_rows = T, cluster_columns = F, name = "Z-scores", bottom_annotation = anno, column_title = "Samples", row_title = "Genes")
dev.off()