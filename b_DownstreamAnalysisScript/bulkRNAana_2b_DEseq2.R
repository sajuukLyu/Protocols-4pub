
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

# analysis
library(DESeq2)

# for more information, please refer to:
# http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# 2. Load data ------------------------------------------------------------

dataMtx <- readRDS("mid/dataMtx.rds")

colnames(dataMtx)
usedMtx <- dataMtx[, c(
  # c("wanted sample names"),
  c(3, 6, 4, 5, 1, 2),
  NULL
)]

# 3. Analyze --------------------------------------------------------------

se <- SummarizedExperiment(assays = list(counts = as.matrix(usedMtx)))
se$condition <- factor(
  rep(c("Fib", "CiPS", "ES"), each = 2),
  levels = c("Fib", "CiPS", "ES")
)

dds <- DESeqDataSet(se, design = ~ condition)

# filter genes with low expression level
rs <- rowMeans(usedMtx)
geneKeep <- rs > quantile(rs, 0.4)
sum(geneKeep)
dds <- dds[geneKeep, ]

# save vsd profile for visualization later
vsd <- vst(dds, blind = T)
saveRDS(vsd, "mid/vsd.rds")

# perform DEG test
dds <- DESeq(dds)

resultsNames(dds)

DEGres <- list()
DEGres$CiPS_vs_Fib <- lfcShrink(dds, coef = "condition_CiPS_vs_Fib", type = "apeglm")
DEGres$ES_vs_Fib <- lfcShrink(dds, coef = "condition_ES_vs_Fib", type = "apeglm")

saveRDS(DEGres, "mid/DEGres.rds")
iwalk(DEGres, ~ write.csv(.x, glue("mid/{.y}.DEG.csv")))

saveRDS(dds, "mid/dds.rds")
