
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
  trimmomaticDir = "app/Trimmomatic-0.39",
  trimAdapter = "TruSeq3-PE-2.fa",
  mapRef = "ref/hg19-STAR-cr3",
  inferDir = "",
  inferGene = "",
  countAnno = "genes.gtf"
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
file_name <- list.files("2_trim", ".gz") # after qc

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

# * * 2.2. Trim -----------------------------------------------------------

dir.create("2_trim")
dir.create(".2_untrim")

trim <- settingList$trimmomaticDir
adapter <- settingList$trimAdapter

trim_cmd <- glue(
  "java -jar {trim}/trimmomatic-0.39.jar PE -threads 10 \\
  0_fastq/{R1} 0_fastq/{R2} \\
  2_trim/{sample_name}_1.trim.fastq.gz \\
  .2_untrim/{sample_name}_1.untrim.fastq.gz \\
  2_trim/{sample_name}_2.trim.fastq.gz \\
  .2_untrim/{sample_name}_2.untrim.fastq.gz \\
  ILLUMINACLIP:{trim}/adapters/{adapter}:2:30:7:1:true \\
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > 2_trim/{sample_name}.trim.log 2>&1")
cat(trim_cmd[1])

trim_dir <- "2_trim"
trim_dir %>% str_c("code/", .) %>% dir.create()

setCMD(trim_cmd, str_c("code/", trim_dir), 1, F)
# multiqc -o 2_trim -f -n trim 2_trim/*.log

# * * 2.3. QC for trimmed data --------------------------------------------

# reload trimmed data first

qc_dir <- "3_qc"
qc_dir %>% dir.create()

qc_cmd <- glue("fastqc -o {qc_dir} -t 8 {trim_dir}/{file_name} &")
cat(qc_cmd[1])

write.table(c("#!/bin/bash\n", qc_cmd), glue("code/{qc_dir}.sh"), quote = F, row.names = F, col.names = F)
# multiqc -o 3_qc -f -n qc.trim 3_qc/*.zip

# * * 2.4. Map ------------------------------------------------------------

ref <- settingList$mapRef

map_dir <- "4_map"
map_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

map_cmd <- glue(
  "STAR --genomeDir {ref} \\
  --runThreadN 20 --outFileNamePrefix {map_dir}/{sample_name} \\
  --readFilesIn {trim_dir}/{R1} {trim_dir}/{R2} \\
  --readFilesCommand zcat --outSAMtype BAM Unsorted \\
  --outSAMstrandField intronMotif \\
  --outFilterIntronMotifs RemoveNoncanonical")
cat(map_cmd[1])

setCMD(map_cmd, str_c("code/", map_dir), 1, T)
# multiqc -o 4_map -f -n map 4_map/*.out

# * * 2.5. Infer ----------------------------------------------------------

infer_path <- settingList$inferDir
gene_bed <- settingList$inferGene

infer_dir <- "5_infer"
infer_dir %>% dir.create()

infer_cmd <- glue(
  "{infer_path}/infer_experiment.py \\
  -r {gene_bed} -i {map_dir}/{sample_name}Aligned.out.bam > {infer_dir}/{sample_name}.infer.txt 2>&1 &")
cat(infer_cmd[1])

write.table(c("#!/bin/bash\n", infer_cmd), glue("code/{infer_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.6. Count ----------------------------------------------------------

gtf <- settingList$countAnno

count_dir <- "6_count"
count_dir %>% dir.create()

count_cmd <- glue(
  "featureCounts -T 20 -g gene_id --extraAttributes gene_name -a {gtf} \\
  -o {count_dir}/counts.txt {samples}> {count_dir}/count.log 2>&1 &",
  samples = str_c(map_dir, "/", sample_name, "Aligned.out.bam \\\n", collapse = ""))
cat(count_cmd)

write.table(c("#!/bin/bash\n", count_cmd), glue("code/{count_dir}.sh"), quote = F, row.names = F, col.names = F)
# multiqc -o 6_count -f -n count 6_count

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

