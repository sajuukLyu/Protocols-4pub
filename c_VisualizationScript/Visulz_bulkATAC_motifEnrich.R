
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
library(ggplot2)
library(ggrepel)
library(scales)
library(latex2exp)
library(patchwork)

# * 2. Load data ----------------------------------------------------------

motifFile <- list.files("motif", "known", full.names = T)
motifGroup <- str_remove_all(motifFile, ".*/|_.*")
motifData <- map(motifFile, fread) %>% set_names(motifGroup)

motifMeta <- fread("../../data/motifTable.txt", sep = "\t")
motifMeta[Symbol == "", Symbol := `Factor Name`][]

name2symbol <- structure(
  motifMeta$Symbol,
  names = motifMeta$Name
)

motifData <- map(motifData, ~ {.x[, symbol := name2symbol[`Motif Name`]]})
motifData <- map(motifData, ~ {.x[is.na(symbol), symbol := str_remove(`Motif Name`, "/.*") %>% str_remove("\\(.*")]})
motifData <- map(motifData, ~ {.x[, symbol := toupper(symbol)]})
motifData <- map(motifData, ~ {.x[, p := -log10(`P-value`)]})

saveRDS(motifData$g1[, .(`Motif Name`, symbol)], "motifAnno.rds")

# * 3. Plot ---------------------------------------------------------------

plotList <- imap(motifData, ~ (
  ggplot(.x, aes(x = seq_along(p), y = p)) +
    geom_point(aes(size = p), show.legend = F) +
    geom_text_repel(
      data = .x[1:15], aes(label = symbol),
      nudge_x = 10, size = 3, fontface = 3, segment.size = .3, segment.alpha = .5) +
    scale_radius(range = c(.1, 2)) +
    labs(x = "order", y = TeX("$-log_{10}(\\textit{P}\\,value)$"), title = .y) +
    theme(
      aspect.ratio = 1,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line()
    )
))

reduce(plotList, `+`) + plot_layout(ncol = 1)

ggsave("graphics/motifEnrich.png", width = 4, height = length(motifData) * 4)

