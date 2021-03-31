
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
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)

# * 2. Load data ----------------------------------------------------------

diffData <- fread("DESeq2/XF_vs_F.DEG.csv")
colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, baseMean := log2(baseMean)][]

# * 3. Plot ---------------------------------------------------------------

dir.create("graphics")

diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene <- diffData[order(p, decreasing = T)][type == "up"][1:10]

pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(diffData, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  labs(
    x = TeX("$log_{2}(base\\,Mean)$"),
    y = TeX("$log_{2}(Fold\\,Change)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())

ggsave("graphics/XF_vs_F.MAplot.png", width = 4, height = 4)
