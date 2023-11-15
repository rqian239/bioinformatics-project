# From https://www.biostars.org/p/350710/

# Use topGO to perform GO enrichment analysis

library(topGO)
library(org.Hs.eg.db)
library(ggplot2)

# Get the significantly differentially expressed genes
diff_expressed_genes <- readr::read_tsv("./results/diff_expression_results.tsv")

# Get gene names
genes <- setNames(diff_expressed_genes$padj, diff_expressed_genes$Gene)

# Get GO annotations
selection <- function(allScore){ return(allScore < 0.05)}
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="SYMBOL")
GOdata <- new("topGOdata",
  ontology="BP",
  allGenes=genes,
  annot=annFUN.GO2genes,
  GO2genes=allGO2genes,
  geneSel=selection,
  nodeSize=10)

# Run Kolmogorov-Smirnov test (as provided by https://www.biostars.org/p/350710/)
results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

# Plot the results
ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
    stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
    xlab("Biological process") +
    ylab("Enrichment") +
    ggtitle("Title") +
    scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
    theme_bw(base_size=24) +
    theme(
        legend.position='none',
        legend.background=element_rect(),
        plot.title=element_text(angle=0, size=24, face="bold", vjust=1),
        axis.text.x=element_text(angle=0, size=18, face="bold", hjust=1.10),
        axis.text.y=element_text(angle=0, size=18, face="bold", vjust=0.5),
        axis.title=element_text(size=24, face="bold"),
        legend.key=element_blank(),     #removes the border
        legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
        legend.text=element_text(size=18),  #Text size
        title=element_text(size=18)) +
    guides(colour=guide_legend(override.aes=list(size=2.5))) +
    coord_flip()



