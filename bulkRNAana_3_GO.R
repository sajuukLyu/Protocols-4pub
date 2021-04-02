
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
#
# require: R whatever
#
# ---

# * 1. Load packages ------------------------------------------------------

setwd("exampleData/RNA")

# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)

# * 2. Load data ----------------------------------------------------------

anno <- readRDS("../../data/hg19anno.rds")

diffData <- fread("DESeq2/XF_vs_F.DEG.csv")
colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]

# * 3. Analyze ------------------------------------------------------------

diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

geneList <- list()
geneList$XFup <- diffData[type == "up", gene]
geneList$Fup <- diffData[type == "down", gene]

egoList <- map(geneList, ~ {
  enrichGO(
    gene = na.omit(select(org.Hs.eg.db, keys = .x, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID),
    OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, readable = T)
})
names(egoList)

iwalk(egoList, ~ write.csv(.x@result, str_c(.y, ".GO.csv")))


