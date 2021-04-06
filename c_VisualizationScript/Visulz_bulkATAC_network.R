
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

# for graph
library(XML)

# graphics
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(scales)
library(ggnewscale)

# * 2. Load data ----------------------------------------------------------

eData <- fread("network/4_TFgene.csv")
vData <- fread("network/4_node.csv")

setkey(vData, Id)

graphData <- xmlParse("network/g4.gexf") %>% xmlToList()

vData[map(graphData$graph$nodes, ".attrs") %>% map("id") %>% unlist() %>% unname(), `:=` (
  x = map(graphData$graph$nodes, "position") %>% map("x") %>% unlist() %>% unname() %>% as.numeric(),
  y = map(graphData$graph$nodes, "position") %>% map("y") %>% unlist() %>% unname() %>% as.numeric()
)][]

eData[, `:=` (x = vData[Source, x], y = vData[Source, y], xend = vData[Target, x], yend = vData[Target, y])][]

# * 3. Plot ---------------------------------------------------------------

vData[, color := "other"][type == "tf", color := cut(Weight, 5)][]

colorName <- vData[type == "tf"][order(-Weight), color] %>% unique() %>% c("other")

geneColor <- structure(
  c(viridis(5, begin = .2, direction = -1), "gray80"),
  names = colorName
)

ggplot() +
  geom_segment(
    data = eData, aes(x, y, xend = xend, yend = yend, size = Weight),
    color = "gray90", alpha = .5, show.legend = F) +
  scale_radius(range = c(0.1, 0.5)) +
  new_scale("size") +
  geom_point(
    data = vData, aes(x, y, color = color, size = Weight),
    show.legend = F) +
  scale_color_manual(values = geneColor) +
  scale_radius(range = c(1, 12)) +
  new_scale("size") +
  geom_text(
    data = vData[color %in% colorName[1:3]], aes(x, y, label = Id, size = Weight),
    fontface = 3, show.legend = F) +
  scale_size(range = c(2, 4)) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank())

ggsave("graphics/g4_network.png", width = 6, height = 6, scale = 1.5)

