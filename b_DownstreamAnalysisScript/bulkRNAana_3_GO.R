
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
#
# require: R whatever
#
# ---

# 1. Load packages --------------------------------------------------------

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

# 2. Load data ------------------------------------------------------------

anno <- readRDS("../../data/hg19anno.rds")

diffList <- readRDS("mid/DEGres.rds")
diffList <- map(diffList, ~ as.data.table(.x, T))
diffList <- map(diffList, ~ {colnames(.x)[1] <- "gene"; .x})
diffList <- map(diffList, ~ .x[is.na(padj), padj := 1])

# 3. Analyze --------------------------------------------------------------

for(i in names(diffList)) {
  message(i)
  
  diffData <- diffList[[i]]
  
  diffData[, type := "ns"]
  diffData[log2FoldChange > 3 & padj < 0.05, type := "up"][log2FoldChange < -3 & padj < 0.05, type := "down"]
  table(diffData$type)
  
  geneList <- list(
    up = diffData[type == "up", gene],
    down = diffData[type == "down", gene]
  )
  
  egoList <- map(geneList, ~ {
    enrichGO(
      gene = na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = .x, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID),
      OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, readable = T)
  })
  
  iwalk(egoList, ~ write.csv(.x@result, str_c("mid/", i, ".", .y, ".GO.csv")))
}


