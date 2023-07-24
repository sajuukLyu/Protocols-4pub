
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

library(tidyverse)
library(magrittr)
library(data.table)

# 2. Load data ------------------------------------------------------------

# This file is the first 7 columns of the featureCount output file.
anno <- readRDS("../../data/hg19anno.rds")

dataPath <- "data"
usedData <- list.files(dataPath, "txt", full = T)

fread(usedData)

# the first 7 cols can be saved as an annotation file.
dataMtx <- map(usedData, ~ fread(.x)[, -1:-7]) %>% reduce(cbind) %>% as.matrix()

rownames(dataMtx) <- uniquifyName(anno$gene_name, anno$Geneid)
colnames(dataMtx) %<>% str_replace_all(".*map/|Aligned.*", "")

saveRDS(dataMtx, "mid/dataMtx.rds")

# Function ----------------------------------------------------------------

uniquifyName <- function(name, ID, sep = "-") {
  
  if(any(duplicated(ID))) {
    stop("Require unique IDs.", call. = F)
  }
  
  dupName <- unique(name[duplicated(name)])
  
  out <- name
  out[name %in% dupName] <- paste0(out[name %in% dupName], sep, ID[name %in% dupName])
  
  out
}

