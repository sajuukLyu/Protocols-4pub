
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
library(rtracklayer)
library(ggplot2)
library(ggsci)
library(patchwork)

# * 2. Load data ----------------------------------------------------------

transGene <- readRDS("../../data/hg19canonicalTrans.rds") %>% as.data.table()
transMeta <- readRDS("../../data/hg19trans.rds") %>% as.data.table()
exonMeta <- readRDS("../../data/hg19exon.rds") %>% as.data.table()

bwFile <- list.files("/mnt/f/exampleData/ATAC/bw", full.names = T)[c(2, 1, 3)]
sampleName <- str_remove_all(bwFile, ".*/|.bw")

# * 3. Plot ---------------------------------------------------------------

plotTrack <- function(gene = "SP5") {
  stopifnot(gene %in% transGene$gene_name)
  
  trans <- transGene[gene_name == gene, trans_id]
  chr <- transMeta[tx_name == trans, seqnames] %>% as.vector()
  startPos <- transMeta[tx_name == trans, start]
  endPos <- transMeta[tx_name == trans, end]
  strand <- transMeta[tx_name == trans, strand] %>% as.vector()
  
  extend <- max(2000, (endPos - startPos) * 0.05)
  startEx <- startPos - extend
  endEx <- endPos + extend
  
  exonData <- exonMeta[name == gene]
  exonData <- exonData[map_lgl(tx_name, ~ {trans %in% .x}), .(start, end)]
  
  usedData <- map(bwFile, ~ import.bw(
    .x, format = "bw", selection = BigWigSelection(GRanges(chr, IRanges(startEx, endEx)))))
  names(usedData) <- sampleName
  
  usedData %<>% map(~ as.data.table(.x))
  usedData %<>% imap(~ {data.table(sample = .y, pos = .x$start[1]:.x$end[nrow(.x)], value = rep(.x$score, .x$width))})
  
  plotData <- purrr::reduce(usedData, rbind)
  plotData$sample %<>% factor(levels = sampleName)
  
  usedColor <- ArchR::ArchRPalettes$stallion %>% unname()
  
  g1 <- ggplot(plotData, aes(pos, value)) +
    geom_area(aes(fill = sample), show.legend = F) +
    scale_fill_manual(values = usedColor) +
    labs(y = "Signal intensity") +
    scale_y_continuous(expand = expansion(c(0, 0.05))) +
    scale_x_continuous(expand = c(0, 0)) +
    facet_wrap(~ sample, ncol = 1, strip.position = "right") +
    theme(
      aspect.ratio = 1/10,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA),
      panel.spacing = unit(0, "lines"),
      strip.background = element_blank(),
      plot.margin = margin(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  segmentData <- data.table(
    start = seq(startPos, endPos, length.out = 10)[2:9]
  )
  
  if(strand == "-") {
    segmentData[, end := start - 1][]
  }else{
    segmentData[, end := start + 1][]
  }
  
  sizeBar <- log10(endPos - startPos) %>% floor() %>% 10^.
  
  g2 <- ggplot(exonData) +
    geom_linerange(x = 0, ymin = startPos, ymax = endPos, size = 0.5) +
    geom_linerange(aes(x = 0, ymin = start, ymax = end), size = 3) +
    geom_segment(
      data = segmentData, aes(y = start, yend = end), x = 0, xend = 0,
      arrow = arrow(length = unit(.4, "lines")), size = 0.3) +
    annotate("segment", x = 1, xend = 1, y = endEx - sizeBar, yend = endEx, size = 1.5) +
    annotate("text", x = 1.2, y = endEx, label = humanFormat(sizeBar), hjust = 1, size = 3) +
    scale_y_continuous(expand = c(0, 0), limits = c(startEx, endEx)) +
    scale_x_continuous(expand = expansion(c(.1, .1))) +
    labs(x = "", y = gene) +
    coord_flip() +
    theme(
      aspect.ratio = 1/10,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.margin = margin()
    )
  
  g1 + g2 + plot_layout(ncol = 1, heights = c(length(sampleName), 1))
  
  ggsave(str_c("graphics/", gene, "_track.png"), width = 10, height = length(sampleName) + 1)
}

humanFormat <- function(n) {
  stopifnot(is.numeric(n), n > 0)
  
  size <- log10(n)
  
  if(size < 3) {
    as.character(n)
  } else if(size < 6) {
    paste0(round(n / 1e3), "K")
  } else if(size < 9) {
    paste0(round(n / 1e6), "M")
  } else {
    paste0(round(n / 1e9), "G")
  }
}

plotTrack("SP5")
plotTrack("CAT")

