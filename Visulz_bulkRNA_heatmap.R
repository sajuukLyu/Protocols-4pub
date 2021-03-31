
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

# graphics
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggsci)
library(scales)

# * 2. Load data ----------------------------------------------------------

vsd <- readRDS("DESeq2/vsd.rds")
vsdMtx <- assay(vsd)

diffData <- fread("DESeq2/XF_vs_F.DEG.csv")
colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]

# * 3. Plot ---------------------------------------------------------------

dir.create("graphics")

diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

upGene <- diffData[order(p, decreasing = T)][type == "up"][1:10, gene]
downGene <- diffData[order(p, decreasing = T)][type == "down"][1:10, gene]

# heatmap with many genes

diffGene <- diffData[type != "ns", gene]
geneType <- diffData[type != "ns", type]
markGene <- c(upGene, downGene)

heatData <- vsdMtx[diffGene, ]
heatData %<>% apply(1, scale) %>% t() %>% set_colnames(colnames(vsdMtx))

colorFun <- colorRamp2(seq(2, -2, len = 9), brewer.pal(9, "RdBu"))

png("graphics/heatmapMany.png", width = 6, height = 8, units = "in", res = 300)

Heatmap(
  heatData, col = colorFun, border = F,
  cluster_rows = T, cluster_columns = T,
  show_row_names = F, show_column_names = T,
  show_row_dend = F, show_column_dend = T,
  column_names_rot = 45,
  row_split = geneType,
  left_annotation = rowAnnotation(
    type = anno_block(
      width = unit(.1, "in"), gp = gpar(fill = pal_nejm()(2)[2:1]))),
  right_annotation = rowAnnotation(
    gene = anno_mark(
      match(markGene, diffGene), markGene,
      labels_gp = gpar(fontface = "italic", fontsize = 10))),
  heatmap_legend_param = list(
    title = "Scaled expression",
    title_position = "lefttop-rot",
    legend_height = unit(2, "in")),
  width = unit(3, "in"),
  height = unit(6, "in"))

dev.off()


# heatmap with few genes

setkey(diffData, gene)

usedGene <- c(upGene, downGene)
geneType <- diffData[usedGene, type]

heatData <- vsdMtx[usedGene, ]
heatData %<>% apply(1, scale) %>% t() %>% set_colnames(colnames(vsdMtx))

colorFun <- colorRamp2(seq(2, -2, len = 9), brewer.pal(9, "RdBu"))

png("graphics/heatmapFew.png", width = 6, height = 6, units = "in", res = 300)

Heatmap(
  heatData, col = colorFun, border = F,
  cluster_rows = T, cluster_columns = T,
  show_row_names = T, show_column_names = T,
  show_row_dend = F, show_column_dend = T,
  column_names_rot = 45,
  row_names_gp = gpar(fontface = "italic"),
  row_split = geneType,
  cell_fun = function(j, i, x, y, width, height, fill){
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "gray40", fill = NA))},
  left_annotation = rowAnnotation(
    type = anno_block(
      width = unit(.1, "in"), gp = gpar(fill = pal_nejm()(2)[2:1]))),
  heatmap_legend_param = list(
    title = "Scaled expression",
    title_position = "lefttop-rot",
    legend_height = unit(2, "in")),
  width = unit(3, "in"),
  height = unit(4, "in"))

dev.off()

