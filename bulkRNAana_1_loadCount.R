
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

library(tidyverse)
library(magrittr)

# * 2. Load data ----------------------------------------------------------

# This file is the first 7 columns of the featureCount output file.
anno <- read_rds("hg19anno.rds")

dataPath <- "data"
usedData <- list.files(dataPath, full = T)

dataMtx <- map(usedData, ~ read_delim(.x, delim = "\t", comment = "#")) %>%
  map(~ .x[-1:-7]) %>% purrr::reduce(cbind) # the first 7 cols can be saved as an annotation file.

rownames(dataMtx) <- scater::uniquifyFeatureNames(anno$Geneid, anno$gene_name)
colnames(dataMtx) %<>% str_replace_all(".*map/|Aligned.*", "")

saveRDS(dataMtx, "dataMtx.rds")

