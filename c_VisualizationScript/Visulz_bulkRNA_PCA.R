
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
library(irlba)

# graphics
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)

# * 2. Load data ----------------------------------------------------------

vsd <- readRDS("DESeq2/vsd.rds")
vsdMtx <- assay(vsd)

diffData <- fread("DESeq2/XF_vs_F.DEG.csv")
colnames(diffData)[1] <- "gene"

# * 3. Plot ---------------------------------------------------------------

dir.create("graphics")

diffData[, type := "ns"][]
diffData[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]
usedGene <- diffData[type %in% c("up", "down"), gene]

PCAdata <- prcomp_irlba(vsdMtx[usedGene, ], n = 3, scale. = T)
PCprop <- (PCAdata$sdev)^2 %>% {round(. / sum(.), 3) * 100}

plotData <- as.data.table(PCAdata$rotation)
plotData[, id := colnames(vsdMtx)][, type := rep(c("F", "XF"), c(3, 4))][]

ggplot(plotData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = type), show.legend = F) +
  geom_text_repel(aes(label = id)) +
  scale_color_d3() +
  labs(
    x = str_c("PC1 (", PCprop[1], "%)"),
    y = str_c("PC2 (", PCprop[2], "%)")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA)
  )

ggsave("graphics/PCA.png", width = 4, height = 4)
