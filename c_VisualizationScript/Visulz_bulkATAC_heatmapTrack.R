
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

# graphics
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggsci)
library(scales)

# * 2. Load data ----------------------------------------------------------

peakData <- fread("peak/peak_value.txt", skip = 3)
peakGroup <- readRDS("peakGroup.rds")

sample <- rep(c("F", "P", "XF"), each = 300) %>% factor(levels = c("P", "F", "XF"))
peakType <- rep(c("g1", "g4"), c(length(peakGroup$`1`), length(peakGroup$`4`))) %>% factor(levels = c("g4", "g1"))

heatData <- as.matrix(peakData)
heatData[is.nan(heatData)] <- 0
max(heatData)

colorFun <- colorRamp2(seq(6, 0, len = 11), brewer.pal(11, "Spectral"))

p <- Heatmap(
  heatData, border = F, col = colorFun,
  cluster_rows = F, cluster_columns = F,
  show_row_names = F, show_column_names = F,
  column_title_side = "bottom",
  row_split = peakType,
  column_split = sample,
  left_annotation = rowAnnotation(
    type = anno_block(
      gp = gpar(col = NA, fill = pal_nejm()(2)[2:1]),
      width = unit(.1, "in"))
  ),
  heatmap_legend_param = list(
    title = "Normalized signal",
    title_position = "lefttop-rot",
    legend_height = unit(2, "in")),
  width = unit(4, "in"),
  height = unit(8, "in")
)

png("graphics/heatmapTrack.png", res = 300, width = 6, height = 9, units = "in")
p
dev.off()

