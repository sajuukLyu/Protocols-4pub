
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
#
# require: R whatever
#
# ---

# * 1. Load packages ------------------------------------------------------

setwd("project/path")

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
usedMtx <- dataMtx[c(
  c("wanted sample names"),
  NULL
)]

# * 3. Analyze ------------------------------------------------------------

dir.create("DESeq2")

se <- SummarizedExperiment(assays = list(counts = as.matrix(usedMtx)))
se$condition <- factor(
  c("a", "b", "c"),
  levels = c("a", "b", "c")
)

dds <- DESeqDataSet(se, design = ~condition)

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

resultsNames(earlyDDS)

DEGres <- list()
DEGres$a_vs_b <- lfcShrink(dds, coef = "condition_a_vs_b", type = "apeglm")
DEGres$a_vs_c <- lfcShrink(dds, coef = "condition_a_vs_c", type = "apeglm")

iwalk(DEGres, ~ write_csv(.x, glue("DESeq2/{.y}.DEG.csv")))
