
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
library(ggplot2)
library(ggchicklet)
library(ggsci)
library(scales)

# * 2. Load data ----------------------------------------------------------

peakAnnoList <- readRDS("peakAnnoList.rds")

# * 3. Plot ---------------------------------------------------------------

plotData <- imap(peakAnnoList, ~ {as.data.table(.x@anno)[, group := .y][, .(annotation, group)]}) %>% purrr::reduce(rbind)

plotData[str_detect(annotation, "Exon"), annotation := "Exon"]
plotData[str_detect(annotation, "Intron"), annotation := "Intron"]
plotData[str_detect(annotation, "Downstream|Intergenic"), annotation := "Intergenic"]
plotData <- plotData[, .N, by = .(annotation, group)]

ggplot(plotData, aes(x = group, y = N)) +
  geom_chicklet(aes(fill = fct_rev(annotation)), position = "fill", width = .95) +
  scale_fill_nejm() +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = "", y = "") +
  theme(
    aspect.ratio = .5,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.title = element_blank()
  )

ggsave("graphics/peakAnnoDisp.png", width = 5, height = 2.5)
