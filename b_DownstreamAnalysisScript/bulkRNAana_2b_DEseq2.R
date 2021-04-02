
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

# analysis
library(DESeq2)

# for more information, please refer to:
# http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# * 2. Load data ----------------------------------------------------------

dataMtx <- readRDS("dataMtx.rds")

colnames(dataMtx)
usedMtx <- dataMtx[, c(
  # c("wanted sample names"),
  1:7,
  NULL
)]

# * 3. Analyze ------------------------------------------------------------

dir.create("DESeq2")

se <- SummarizedExperiment(assays = list(counts = as.matrix(usedMtx)))
se$condition <- factor(
  rep(c("F", "XF"), c(3, 4)),
  levels = c("F", "XF")
)

dds <- DESeqDataSet(se, design = ~ condition)

# filter genes with low expression level
rs <- rowMeans(usedMtx)
geneKeep <- rs > quantile(rs, 0.4)
sum(geneKeep)
dds <- dds[geneKeep, ]

# save vsd profile for visualization later
vsd <- vst(dds, blind = T)
saveRDS(vsd, "DESeq2/vsd.rds")

# perform DEG test
dds <- DESeq(dds)

resultsNames(dds)

DEGres <- list()
DEGres$XF_vs_F <- lfcShrink(dds, coef = "condition_XF_vs_F", type = "apeglm")

iwalk(DEGres, ~ write.csv(.x, glue("DESeq2/{.y}.DEG.csv")))

saveRDS(dds, "DESeq2/dds.rds")
