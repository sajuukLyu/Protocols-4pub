
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
#
# require: R whatever
#
# ---

# * 1. Load packages ------------------------------------------------------

setwd("exampleData/ATAC")

# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(org.Hs.eg.db)
library(clusterProfiler)

# * 2. Load data ----------------------------------------------------------

anno <- readRDS("../../data/hg19anno.rds")
peakAnnoList <- readRDS("peakAnnoList.rds")

id2gene <- structure(anno$gene_name, names = anno$Geneid)

annoGene <- map(peakAnnoList, ~ {.x@anno$geneId})

egoList <- map(annoGene, ~ {
  enrichGO(
    gene = na.omit(select(org.Hs.eg.db, keys = .x, columns = "ENTREZID", keytype = "ENSEMBL")$ENTREZID),
    OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, readable = T)
})
names(egoList)

iwalk(egoList, ~ write.csv(.x@result, str_c("peakGroup_", .y, ".GO.csv")))

