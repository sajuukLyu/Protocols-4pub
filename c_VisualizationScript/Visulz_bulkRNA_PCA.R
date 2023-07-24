
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
library(irlba)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)

# 2. Load data ------------------------------------------------------------

vsd <- readRDS("mid/vsd.rds")
vsdMtx <- assay(vsd)

diffData <- fread("mid/ES_vs_Fib.DEG.csv")
colnames(diffData)[1] <- "gene"

# 3. Plot -----------------------------------------------------------------

diffData[, type := "ns"][]
diffData[is.na(padj), padj := 1][]
diffData[log2FoldChange > 3 & padj < 0.05, type := "up"][log2FoldChange < -3 & padj < 0.05, type := "down"][]
usedGene <- diffData[type %in% c("up", "down"), gene]

PCAdata <- prcomp_irlba(t(vsdMtx[usedGene, ]), n = 3, scale. = T)

s <- summary(PCAdata)
s <- s$importance[2, ] %>% round(4)

plotData <- as.data.table(PCAdata$x)
plotData[, id := colnames(vsdMtx)][, type := rep(c("Fib", "CiPS", "ES"), each = 2)][]

usedCol <- pal_npg()(10)[c(1, 4, 3)] %>% set_names(c("Fib", "CiPS", "ES"))

ggplot(plotData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = type), size = 3, shape = 16, show.legend = F) +
  geom_text_repel(aes(label = id)) +
  scale_color_manual(values = usedCol) +
  labs(
    x = str_c("PC1 (", s[1] * 100, "%)"),
    y = str_c("PC2 (", s[2] * 100, "%)")) +
  coord_fixed() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA)
  )

ggsave("plot/PCA.png", width = 6, height = 2)
