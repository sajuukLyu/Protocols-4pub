
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

library(tidyverse)
library(magrittr)

# * 2. Load data ----------------------------------------------------------

# This file is the first 7 columns of the featureCount output file.
anno <- readRDS("../../data/hg19anno.rds")

dataPath <- "data"
usedData <- list.files(dataPath, "txt", full = T)

dataMtx <- map(usedData, ~ read_delim(.x, delim = "\t", comment = "#")) %>%
  map(~ .x[-1:-7]) %>%
  purrr::reduce(cbind) %>%
  as.matrix() # the first 7 cols can be saved as an annotation file.

rownames(dataMtx) <- scater::uniquifyFeatureNames(anno$Geneid, anno$gene_name)
colnames(dataMtx) %<>% str_replace_all(".*map/|Aligned.*", "")

saveRDS(dataMtx, "dataMtx.rds")

