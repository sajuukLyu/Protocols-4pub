
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
library(ChIPseeker)
library(GenomicFeatures)

# * 2. Load data ----------------------------------------------------------

txdb <- makeTxDbFromGFF("../../data/hg19genes.gtf")

peakMeta <- readRDS("peakMeta.rds")
peakMeta[, id := paste(Chr, Start, End, sep = "_")][]

peakIndex <- map(peakGroup, ~ {match(.x, peakMeta$id)})

# * 3. Analyze ------------------------------------------------------------

peakGRlist <- map(peakIndex, ~ {
  GRanges(
    seqnames = peakMeta[.x, Chr],
    ranges = IRanges(
      start = peakMeta[.x, Start],
      end = peakMeta[.x, End]))})

peakAnnoList <- map(peakGRlist, annotatePeak, TxDb = txdb, tssRegion = c(-2000, 2000))

map(peakGRlist, ~ {as.data.table(.x)[, c("width", "strand") := NULL][, id := 1:.N][, seqnames := str_c("chr", seqnames)]}) %>%
  iwalk(~ fwrite(.x, str_c("peakGroup_", .y, ".bed"), sep = "\t", col.names = F))

saveRDS(peakAnnoList, "peakAnnoList.rds")
