
# MESSAGE -----------------------------------------------------------------
# 
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
# 
# require: R whatever
# 
# ---

# * 0. Setting ------------------------------------------------------------

settingList <- list(
  workDir = ".",
  mapRef = "fasta",
  lambdaRef = "ref/lambda",
  BSfilter = "BSfilter.R",
  bdg2bwDir = "app/bedGraphToBigWig",
  refLen = "fasta/genome.fa.fai"
)

# * 1. Load packages ------------------------------------------------------

setwd(settingList$workDir)

library(tidyverse)
library(magrittr)
library(glue)

# * 2. Preprocess ---------------------------------------------------------

setwd("0_fastq")

list.files(".")

from_file <- list.files(".", "")
to_file <- gsub("", "", from_file)

file.rename(from_file, to_file)

setwd("..")

# * * 2.0. Load data ------------------------------------------------------

file_name <- list.files("0_fastq", ".gz")

R1 <- grep("_1\\.", file_name, value = T)
R2 <- grep("_2\\.", file_name, value = T)

sample_name <- gsub("_1\\..*", "", R1)

# * * 2.1. QC for raw data ------------------------------------------------

dir.create("code")

qc_dir <- "1_qc"
qc_dir %>% dir.create()

qc_cmd <- glue("fastqc -o {qc_dir} -t 8 0_fastq/{file_name} &")
cat(qc_cmd[1])

write.table(c("#!/bin/bash\n", qc_cmd), glue("code/{qc_dir}.sh"), quote = F, row.names = F, col.names = F)
# multiqc -o 1_qc -f -n qc.raw 1_qc/*.zip

# * * 2.2. Map ------------------------------------------------------------

ref <- settingList$mapRef

file_dir <- "0_fastq"
map_dir <- "2_map"
map_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

map_cmd <- glue(
  "bismark --multicore 5 -X 700 --dovetail -p 2 --un --ambiguous --genome {ref} \\
  {file_dir}/{R1} {file_dir}/{R2} \\
  {map_dir} > {map_dir}/{sample_name}.log")
cat(map_cmd[1])

setCMD(map_cmd, str_c("code/", map_dir), 1, T)
# multiqc -o 2_map -f -n map 2_map/*.txt

# * * 2.3. Dedup ----------------------------------------------------------

dedup_dir <- "3_dedup"
dedup_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

dedup_cmd <- glue(
  "deduplicate_bismark -p -o {sample_name} --output_dir {dedup_dir} \\
  {map_dir}/{sample_name}_1_bismark_bt2_pe.bam")
cat(dedup_cmd[1])

setCMD(dedup_cmd, str_c("code/", dedup_dir), 1, T)

# * * 2.4. Meth -----------------------------------------------------------

meth_dir <- "4_meth"
meth_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

meth_cmd <- glue(
  "bismark_methylation_extractor -p --no_overlap --ignore 5 --ignore_r2 5 \\
  --multicore 7 --bedGraph --cytosine_report --zero_based --CX --buffer_size 50% \\
  --genome_folder {ref} -o {meth_dir} {dedup_dir}/{sample_name}.deduplicated.bam")
cat(meth_cmd[1])

setCMD(meth_cmd, str_c("code/", meth_dir), 1, T)

# * * 2.5. Lambda ---------------------------------------------------------

lambda_ref <- settingList$lambdaRef

lambda_dir <- "5_lambda"
lambda_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

lambda_cmd <- glue(
  "bismark --multicore 5 -X 700 --dovetail -p 2 --un --ambiguous --genome {lambda_ref} \\
  {file_dir}/{R1} {file_dir}/{R2} \\
  {lambda_dir} > {lambda_dir}/{sample_name}.log
  
  deduplicate_bismark -p -o {sample_name} --output_dir {lambda_dir} \\
  {lambda_dir}/{sample_name}_1_bismark_bt2_pe.bam
  
  bismark_methylation_extractor -p --no_overlap --ignore 5 --ignore_r2 5 \\
  --multicore 7 --bedGraph --cytosine_report --zero_based --CX --buffer_size 50% \\
  --genome_folder {lambda_ref} -o {lambda_dir} {lambda_dir}/{sample_name}.deduplicated.bam")

setCMD(lambda_cmd, str_c("code/", lambda_dir), 1, T)

# * * 2.6. Extract --------------------------------------------------------

extract_dir <- "6_extract"
extract_dir %>% dir.create()

extract_cmd <- glue(
  "{awk_cmd} {meth_dir}/{sample_name}.deduplicated.CX_report.txt > {extract_dir}/{sample_name}.cg.txt &",
  awk_cmd = "awk -F '\\t' '{if($4 + $5 >= 5 && $6 == \"CG\"){print $0}}'")
cat(extract_cmd[1])

write.table(c("#!/bin/bash\n", extract_cmd), glue("code/{extract_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.7. Filter ---------------------------------------------------------

p_files <- list.files(lambda_dir, full.names = T)
p_text <- map(p_files, read_delim, delim = "\n")
m <- p_text %>% map(~ .x[[1]][13:15]) %>% map(~ str_replace(.x, ".*\\t", "")) %>% map(~ sum(as.numeric(.x)))
u <- p_text %>% map(~ .x[[1]][16:18]) %>% map(~ str_replace(.x, ".*\\t", "")) %>% map(~ sum(as.numeric(.x)))
p <- map2_dbl(m, u, ~ {.x / (.x + .y)})

filter_dir <- "7_filter"
filter_dir %>% dir.create()

BSfilter <- settingList$BSfilter

filter_cmd <- glue(
  "Rscript {BSfilter} -i {extract_dir}/{sample_name}.cg.txt \\
  -d {filter_dir}/{sample_name}.dss.txt \\
  -b {filter_dir}/{sample_name}.bdg \\
  -p {p} > {filter_dir}/{sample_name}.log 2>&1 &")
cat(filter_cmd[1])

write.table(c("#!/bin/bash\n", filter_cmd), glue("code/{filter_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.8 Bigwig ----------------------------------------------------------

bw_dir <- "8_bw"
bw_dir %>% dir.create()

bdg2bw <- settingList$bdg2bwDir
ref_len <- settingList$refLen

bw_cmd <- glue("{bdg2bw} {filter_dir}/{sample_name}.bdg {ref_len} {bw_dir}/{sample_name}.bw &")
cat(bw_cmd[1])

write.table(c("#!/bin/bash\n", bw_cmd), glue("code/{bw_dir}.sh"), quote = F, row.names = F, col.names = F)

# * 3. Function -----------------------------------------------------------

setCMD <- function(cmd, dir = ".", sepN = 1, clu = F) {
  cmd %>% tapply(seq_along(.) %% sepN, c) %>% imap(~ {
    ifelse(clu, glue(
      "#!/bin/bash
      #SBATCH -J batch{.y}
      #SBATCH -o batch{.y}.%j.out
      #SBATCH -e batch{.y}.%j.err
      #SBATCH -p cn-long
      #SBATCH -N 1
      #SBATCH --ntasks-per-node=20
      #SBATCH --no-requeue
      #SBATCH -A hkdeng_g1
      #SBATCH --qos=hkdengcnl
      #export PATH=/gpfs1/hkdeng_pkuhpc/lvyl/app/Bismark-0.22.3:$PATH
      export PATH=/gpfs1/hkdeng_pkuhpc/lvyl/app/anaconda3/envs/lvyl/bin:$PATH"),
      "#!/bin/bash") %>%
      c(.x)}) %T>%
    iwalk(~ write.table(.x, glue("{dir}/batch{.y}.sh"), quote = F, row.names = F, col.names = F)) %>%
    names() %>% map_chr(~ glue("{head} {dir}/batch{.x}.sh {tail}",
                               head = ifelse(clu, "pkubatch", "sh"),
                               tail = ifelse(clu, "; sleep 1", "&"))) %>%
    c("#!/bin/bash", .) %>% as_tibble() %>%
    write_delim(glue("{dir}/submit.sh"), "\n", col_names = F)
}
