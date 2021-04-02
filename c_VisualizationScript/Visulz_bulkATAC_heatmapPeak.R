
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
library(GenomicRanges)

# graphics
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggsci)
library(scales)

# * 2. Load data ----------------------------------------------------------

diffList <- readRDS("diffList.rds")
peakMeta <- readRDS("peakMeta.rds")
peakMtx <- readRDS("peakMtx.rds")

diffPeak <- map(diffList, ~ {as.data.table(.x)[abs(Fold) > 1, paste(seqnames, start, end, sep = "_")]}) %>%
  unlist() %>% unname() %>% unique()

peakMeta[, id := paste(Chr, Start, End, sep = "_")][]

diffMtx <- as.matrix(peakMtx[peakMeta$id %in% diffPeak]) %>% set_rownames(peakMeta[id %in% diffPeak, id])

# * 3. Plot ---------------------------------------------------------------

dir.create("graphics")

heatData <- apply(diffMtx, 1, scale) %>% t() %>% set_colnames(colnames(diffMtx))

colorFun <- colorRamp2(seq(2, -2, len = 9), brewer.pal(9, "RdBu"))

p <- Heatmap(
  heatData[, c(3:5, 1:2, 6:8)], col = colorFun, border = F,
  cluster_rows = T, cluster_columns = F,
  show_row_names = F, show_column_names = T,
  show_row_dend = F, show_column_dend = F,
  column_names_rot = 45,
  row_km = 5,
  row_km_repeats = 100,
  left_annotation = rowAnnotation(
    type = anno_block(
      width = unit(.1, "in"), gp = gpar(fill = pal_nejm()(5)))),
  heatmap_legend_param = list(
    title = "Scaled expression",
    title_position = "lefttop-rot",
    legend_height = unit(2, "in")),
  width = unit(3, "in"),
  height = unit(6, "in"))

png("graphics/heatmapPeak.png", width = 4, height = 7, units = "in", res = 300)
p
dev.off()

peakOrder <- row_order(p)
peakGroup <- map(peakOrder, ~ {rownames(diffMtx)[.x]})

saveRDS(peakGroup, "peakGroup.rds")
