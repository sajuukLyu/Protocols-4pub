
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
#
# require: R < 4.0 (recommended)
#
# note: DESeq is only for visualization usage in my protocols.
#       It was no longer available since bioconductor 3.11 for R4.
#       The last version of DESeq was built on R3.
#       So this package may not work on R4.
#
# ---

# * 1. Load packages ------------------------------------------------------

setwd("project/path")

# grammar
library(tidyverse)
library(magrittr)
library(glue)

# analysis
library(DESeq)

# for more information, please refer to:
# http://bioconductor.org/packages/3.10/bioc/vignettes/DESeq/inst/doc/DESeq.pdf

# * 2. Load data ----------------------------------------------------------

dataMtx <- readRDS("dataMtx.rds")

colnames(dataMtx)
usedMtx <- dataMtx[c(
  c("wanted sample names"),
  NULL
)]

# * 3. Analyze ------------------------------------------------------------

dir.create("DESeq")

condition <- factor(
  c("a", "b", "c"),
  levels = c("a", "b", "c")
)

cds <- newCountDataSet(usedMtx, condition) %>%
  estimateSizeFactors() %>%
  estimateDispersions(method = "blind", sharingMode = "fit-only")

# save vsd profile for visualization later
vsd <- varianceStabilizingTransformation(cds)
saveRDS(vsd, "DESeq/vsd.rds")

# perform DEG test (BETTER DO NOT, only if there is no other choice)
# DEGres <- list()
# DEGres$a_vs_b <- nbinomTest(cds, "a", "b")
# DEGres$a_vs_c <- nbinomTest(cds, "a", "c")

# iwalk(DEGres, ~ write_csv(.x, glue("DESeq/{.y}.DEG.csv")))
