
# MESSAGE -----------------------------------------------------------------
# 
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
# 
# require: R whatever
# 
# ---

# 0. Load packages --------------------------------------------------------

library(jsonlite)
library(tidyverse)
library(magrittr)
library(glue)

setwd("/mnt/d/proj/github/Protocols-4pub")

settingList <- fromJSON("a_PreprocessPipeline/bulkRNApre_conf_X.json")
cmdConf <- fromJSON("a_PreprocessPipeline/cmd_conf_X.json")

setwd(settingList$workDir)

# A. Rename ---------------------------------------------------------------

setwd("0_fastq")

list.files(".")

from_file <- list.files(".", "")
to_file <- gsub("", "", from_file)

file.rename(from_file, to_file)

setwd("..")

# B. Load data ------------------------------------------------------------

file_name <- list.files("0_fastq", ".gz")
file_name <- list.files("2_trim", ".gz") # after qc

R1 <- grep("_1\\.", file_name, value = T)
R2 <- grep("_2\\.", file_name, value = T)

sample_name <- gsub("_1\\..*", "", R1)

# 1. QC for raw data ------------------------------------------------------

dir.create("code")

qc_dir <- "1_qc"
qc_dir %>% dir.create()
qc_dir %>% str_c("code/", .) %>% dir.create()

qc <- settingList$fastqcDir

qc_cmd <- glue("{qc} -o {qc_dir} -t 48 0_fastq/{file_name}")
cat(qc_cmd[1])

setCMD(qc_cmd, str_c("code/", qc_dir), 6)
# multiqc -o 1_qc -f -n qc.raw 1_qc/*.zip

# 2. Trim -----------------------------------------------------------------

trim_dir <- "2_trim"
trim_dir %>% dir.create()
trim_dir %>% str_c("code/", .) %>% dir.create()
dir.create(".2_untrim")

trim <- settingList$trimmomaticDir
adapter <- settingList$trimAdapter

trim_cmd <- glue(
  "java -jar {trim} PE -threads 48 \\
  0_fastq/{R1} 0_fastq/{R2} \\
  2_trim/{sample_name}_1.trim.fastq.gz \\
  .2_untrim/{sample_name}_1.untrim.fastq.gz \\
  2_trim/{sample_name}_2.trim.fastq.gz \\
  .2_untrim/{sample_name}_2.untrim.fastq.gz \\
  ILLUMINACLIP:{adapter}:2:30:7:1:true \\
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > \\
  2_trim/{sample_name}.trim.log 2>&1"
)
cat(trim_cmd[1])

setCMD(trim_cmd, str_c("code/", trim_dir), 6)
# multiqc -o 2_trim -f -n trim 2_trim/*.log

# 3. QC for trimmed data --------------------------------------------------

# reload trimmed data first, or
file_name <- glue("{rep(sample_name, each = 2)}_{rep(1:2, 6)}.trim.fastq.gz")
R1 <- grep("_1\\.", file_name, value = T)
R2 <- grep("_2\\.", file_name, value = T)

qc_dir <- "3_qc"
qc_dir %>% dir.create()
qc_dir %>% str_c("code/", .) %>% dir.create()

qc <- settingList$fastqcDir

qc_cmd <- glue("{qc} -o {qc_dir} -t 48 {trim_dir}/{file_name}")
cat(qc_cmd[1])

setCMD(qc_cmd, str_c("code/", qc_dir), 6)
# multiqc -o 3_qc -f -n qc.trim 3_qc/*.zip

# 4. Map ------------------------------------------------------------------

map_dir <- "4_map"
map_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

star <- settingList$starDir
ref <- settingList$mapRef

map_cmd <- glue(
  "{star} --genomeDir {ref} \\
  --runThreadN 48 \\
  --outFileNamePrefix {map_dir}/{sample_name} \\
  --readFilesIn {trim_dir}/{R1} {trim_dir}/{R2} \\
  --readFilesCommand zcat \\
  --outSAMtype BAM Unsorted \\
  --outSAMstrandField intronMotif \\
  --outFilterIntronMotifs RemoveNoncanonical"
)
cat(map_cmd[1])

setCMD(map_cmd, str_c("code/", map_dir), 6)
# multiqc -o 4_map -f -n map 4_map/*.out

# 5. Infer (optional) -----------------------------------------------------

infer_dir <- "5_infer"
infer_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

infer_path <- settingList$inferDir
gene_bed <- settingList$inferGene

infer_cmd <- glue(
  "{infer_path} -r {gene_bed} \\
  -i {map_dir}/{sample_name}Aligned.out.bam > \\
  {infer_dir}/{sample_name}.infer.txt 2>&1"
)
cat(infer_cmd[1])

setCMD(infer_cmd, str_c("code/", infer_dir), 6)

# 6. Count ----------------------------------------------------------------

count_dir <- "6_count"
count_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

count_path <- settingList$countDir
gtf <- settingList$countAnno

count_cmd <- glue(
  "{count_path} -s 0 -p -T 48 \\
  -g gene_id --extraAttributes gene_name \\
  -a {gtf} \\
  -o {count_dir}/counts.txt {samples} > \\
  {count_dir}/count.log 2>&1",
  samples = str_c(map_dir, "/", sample_name, "Aligned.out.bam \\\n", collapse = "")
)
cat(count_cmd)

setCMD(count_cmd, str_c("code/", count_dir), 1)
# multiqc -o 6_count -f -n count 6_count

# C. Function -------------------------------------------------------------

setCMD <- function(cmd, dir = ".", sepN = 1, conf = cmdConf) {
  
  idx <- seq_along(cmd) %% sepN
  idx[idx == 0] <- sepN
  
  cmdList <- tapply(cmd, idx, c)
  
  cmdHead <- glue(cmdConf$head, .trim = F)
  cmdList <- map2(cmdHead, cmdList, ~ c(.x, .y))
  names(cmdList) <- 1:length(cmdList)
  
  iwalk(cmdList, ~ write.table(
    .x,
    glue("{dir}/batch{.y}.sh"),
    sep = "\n", quote = F, row.names = F, col.names = F,
    eol = "\n"
  ))
  
  submit <- glue(
    "{cmdConf$prefix} {dir}/batch{names(cmdList)}.sh {cmdConf$suffix}"
  )
  write.table(
    c("#!/bin/bash", submit),
    glue("{dir}/submit.sh"),
    sep = "\n", quote = F, row.names = F, col.names = F,
    eol = "\n"
  )
  
  system(glue("ls -lh {dir}"))
}

