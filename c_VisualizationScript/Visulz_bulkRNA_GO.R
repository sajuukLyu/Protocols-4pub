
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
library(org.Hs.eg.db)
library(clusterProfiler)

# graphics
library(ggplot2)
library(ggsci)
library(latex2exp)
library(patchwork)

# * 2. Load data ----------------------------------------------------------

GOfile <- list.files(".", "GO.csv")

# * 3. Plot ---------------------------------------------------------------

GOdata <- map(GOfile, ~ fread(.x)) %>% set_names(str_remove(GOfile, ".GO.*"))
plotData <- GOdata %>% imap(~ {.x[pvalue < 0.05 & Count >= 5][1:20][, type := .y]})
plotData %<>% map(~ {.x[, p := -log10(pvalue)]})

GOplot <- imap(plotData, ~ {
  ggplot(.x, aes(x = p, y = fct_reorder(Description, p))) +
    geom_col(aes(alpha = p), fill = "black", width = .8, show.legend = F) +
    scale_alpha_continuous(range = c(.5, 1)) +
    scale_x_continuous(expand = expansion(c(0, 0.05))) +
    labs(x = TeX("$-log_{10}(\\textit{P}\\,value)$"), y = "", title = .y) +
    theme(
      aspect.ratio = 0.75,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line())
})

purrr::reduce(GOplot, `+`) + plot_layout(ncol = 1)

ggsave("graphics/GO.png", width = 10, height = 4 * length(GOfile))
