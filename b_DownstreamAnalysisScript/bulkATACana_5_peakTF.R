
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

library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# * 2. Load data ----------------------------------------------------------

peakAnno <- fread("motif/allPeak.txt")
colnames(peakAnno)[1] %<>% str_remove(" .*")

saveRDS(peakAnno[, PeakID:`GC%`], "peakAnnoHomer.rds")
peakGene <- readRDS("peakAnnoHomer.rds")
peakGene <- peakGene[, c("PeakID", "Annotation", "Gene Name")]
colnames(peakGene)[3] <- "gene"
peakGene[, Annotation := str_remove(Annotation, " \\(.*")][]
peakGene[Annotation == "Intergenic", gene := "no"][]
setkey(peakGene, PeakID)

motifAnno <- readRDS("motifAnno.rds")

name2symbol <- structure(
  motifAnno$symbol,
  names = motifAnno$`Motif Name`
)

peakGroup <- readRDS("peakGroup.rds")

# * 3. Analyze ------------------------------------------------------------

selCols <- peakAnno[, Chr:`GC%`] %>% colnames()
peakAnno[, (selCols) := NULL]

motifMeta <- name2symbol[str_remove(colnames(peakAnno)[2:ncol(peakAnno)], " .*")]

motifMeta <- data.table(
  name = names(motifMeta),
  symbol = motifMeta
)

motifMeta[, rep := 1:.N, by = symbol][]
motifMeta[rep != 1, symbol := str_c(symbol, "_", rep)][]
colnames(peakAnno)[2:ncol(peakAnno)] <- motifMeta$symbol

peakTFmtx <- melt(peakAnno, id.vars = "PeakID", variable.name = "TF", value.name = "val")
peakTFbind <- peakTFmtx[val > -log(1e-4)]
peakTFbind[, .N, by = PeakID]$N %>% table()
peakTFbind[, .N, by = PeakID]$N %>% hist(breaks = 50)
setkey(peakTFbind, PeakID)

peakGroup[c("1", "4")] %>% imap(~ {
  n_TFpeak <- peakTFbind[.x][, gene := peakGene[PeakID, gene]][!is.na(val)]
  n_TFgene <- n_TFpeak[!gene %in% c("no", ""), .(val = sum(val)), by = .(TF, gene)] %>% set_colnames(c("Source", "Target", "Weight"))
  n_TFgene[, Weight := rescale(Weight, to = c(0.05, 1))]
  n_node <- table(c(as.vector(n_TFgene$Source), n_TFgene$Target)) %>% as.data.table() %>% set_colnames(c("Id", "Weight"))
  n_node[, `:=` (Label = Id, type = ifelse(Id %in% n_TFgene$Source, "tf", "gene"))]
  
  fwrite(n_TFgene, str_c("network/", .y, "_TFgene.csv"), sep = ",")
  fwrite(n_node, str_c("network/", .y, "_node.csv"), sep = ",")
})

