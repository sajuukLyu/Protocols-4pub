
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
#
# require: R > 4.0
#
# ---

# * 1. Load packages ------------------------------------------------------

setwd("F:/exampleData/ATAC")

# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(DiffBind)

# * 2. Load data ----------------------------------------------------------

bamFile <- list.files("bam", ".bam$", full.names = T)
usedSample <- str_remove_all(bamFile, ".*bam/|.filter.*")
sampleType <- str_remove(usedSample, "-[123]$")
sampleRep <- str_extract(usedSample, "[123]$") %>% as.numeric()

sampleSheet <- data.table(
  SampleID = usedSample,
  Condition = sampleType,
  Replicate = sampleRep,
  bamReads = bamFile,
  Peaks = str_c("peak/", sampleType, ".narrowPeak"),
  PeakCaller = "narrow"
)

dbaAll <- dba(sampleSheet = sampleSheet, minOverlap = 2)

dbaAll$chrmap %<>% str_c("chr", .)
dbaAll <- dba.blacklist(dbaAll, blacklist = DBA_BLACKLIST_HG19, greylist = F)
dbaAll$chrmap %<>% str_remove("chr")

dbaAll <- dba.count(dbaAll, minOverlap = 2, fragmentSize = 200, summits = 250)
dbaAll <- dba.normalize(dbaAll, background = T, normalize = DBA_NORM_NATIVE)
dbaAll <- dba.contrast(dbaAll, minMembers = 2, categories = DBA_CONDITION)
dbaAll <- dba.analyze(dbaAll, bBlacklist = F, bGreylist = F)

saveRDS(dbaAll, "dbaAll.rds")
