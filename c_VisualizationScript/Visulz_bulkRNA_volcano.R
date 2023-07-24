
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

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)

# 2. Load data ------------------------------------------------------------

diffData <- fread("mid/ES_vs_Fib.DEG.csv")
colnames(diffData)[1] <- "gene"

diffData[is.na(padj), padj := 1][]
diffData[, p := -log10(padj)][]

# 3. Plot -----------------------------------------------------------------

diffData[, type := "ns"][]
diffData[log2FoldChange > 3 & padj < 0.05, type := "up"][log2FoldChange < -3 & padj < 0.05, type := "down"][]

labelGene <- diffData[order(p, decreasing = T)][type == "up"][1:10]

pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(diffData, aes(x = log2FoldChange, y = p)) +
  geom_point(aes(color = type, size = p), show.legend = F) +
  geom_hline(yintercept = -log10(0.05), color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 3, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = -3, color = "gray60", linetype = "dashed") +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = "italic",
    nudge_x = .5, nudge_y = .5,
    max.overlaps = 100) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
  labs(
    x = TeX("$log_{2}(Fold\\,Change)$"),
    y = TeX("$-log_{10}(\\textit{P}\\,value)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line()
  )

ggsave("plot/ES_vs_Fib.volcano.pdf", width = 4, height = 4)
