library(topGO)
library(org.Hs.eg.db)
library(grid)
library(gridExtra)
library(insight)

df <- readr::read_tsv("./SRP073813_diff_expressed_genes.tsv")

gene_list <- df$padj
names(gene_list) <- df$Gene

topDiffGenes <- function (allscore) {
    return(allscore > 0.00)
}

sampleGOdata <- new("topGOdata", description = "simple", ontology = "BP",
                        allGenes = gene_list, geneSelectionFun = topDiffGenes,
                        annot = annFUN.org, mapping = "org.Hs.eg.db", 
                        ID = "symbol")

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")

elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                    classicKS = resultKS, elimKS = elim,
                    orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

grid.table(allRes)

pdf('nodes_seq.pdf')
showSigOfNodes(sampleGOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')
dev.off()


pValue.classic <- score(resultKS)
pValue.elim <- score(elim)[names(pValue.classic)]
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

gCol <- colMap(gstat$Significant)

pdf('topgo_seq')
plot(pValue.classic, pValue.elim, xlab = "p-value classic",
        ylab = "p-value elim", pch = 19, cex = gSize, col = gCol)
dev.off()

sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
        elim = pValue.elim[sel.go],
        classic = pValue.classic[sel.go])

export_table(allRes, format="html")
