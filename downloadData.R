
# MESSAGE -----------------------------------------------------------------
# 
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
# 
# require: R whatever
# 
# ---

# * Load packages ---------------------------------------------------------

library(tidyverse)
library(magrittr)
library(glue)

# * From GEO --------------------------------------------------------------

setwd("/mnt/f") # download dir

# RNA
setwd("exampleData/RNA/raw")
geo <- "GSE147839"
sra <- "SRP254790"
bioProj <- "PRJNA622253"

# ATAC
setwd("exampleData/ATAC/raw")
geo <- "GSE157237"
sra <- "SRP279550"
bioProj <- "PRJNA660602"

# download from SRA directly using wget (slow)

dir.create("sra")

download.file(
  glue("http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={sra}"),
  glue("{where}/{sra}.csv", where = getwd()))

projMeta <- read_csv(glue("{where}/{sra}.csv", where = getwd()))

srr <- projMeta$Run
downPath <- projMeta$download_path

down_cmd <- glue("wget -b -c -o sra/{srr}.log -O sra/{srr}.sra {downPath}")

write.table(c("#!/bin/bash\n", down_cmd), "down.sh", sep = "\n", quote = F, row.names = F, col.names = F)

# download from ENA using aspera (fast)
# NOTE: NOT all datasets in SRA are accessible in ENA

# <todo> I will complete this part when needed in the future.


# * From ArrayExpress -----------------------------------------------------

# <todo> I will complete this part when needed in the future.


# * Extract fastq ---------------------------------------------------------

dir.create("fastq")

ext_cmd <- glue("fasterq-dump -e 10 -3 -O fastq sra/{srr}.sra")

write.table(c("#!/bin/bash\n", ext_cmd), "ext.sh", sep = "\n", quote = F, row.names = F, col.names = F)

gzip_cmd <- glue(
  "gzip fastq/{srr}.sra_1.fastq &
  gzip fastq/{srr}.sra_2.fastq &")

write.table(c("#!/bin/bash\n", gzip_cmd), "gzip.sh", sep = "\n", quote = F, row.names = F, col.names = F)

