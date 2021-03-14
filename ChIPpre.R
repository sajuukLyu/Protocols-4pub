
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
library(glue)

# * 2. Preprocess ---------------------------------------------------------

setwd("0_fastq") # contains raw fastq files

list.files(".")

from_file <- list.files(".", "pattern") # rename files if necessary
to_file <- gsub(".*-1a", "sample_names", from_file)

file.rename(from_file, to_file)

setwd("..")

# * * 2.0. Load data ------------------------------------------------------

file_name <- list.files("0_fastq", ".gz")
file_name <- list.files("2_trim", ".gz")

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

trim <- "/lustre/user/liclab/lvyl/app/Trimmomatic-0.39" # path to trimmomatic

trim_cmd <- glue(
  "java -jar {trim}/trimmomatic-0.39.jar PE -threads 10 \\
  0_fastq/{R1} 0_fastq/{R2} \\
  2_trim/{sample_name}_1.trim.fastq.gz \\
  .2_untrim/{sample_name}_1.untrim.fastq.gz \\
  2_trim/{sample_name}_2.trim.fastq.gz \\
  .2_untrim/{sample_name}_2.untrim.fastq.gz \\
  ILLUMINACLIP:{trim}/adapters/TruSeq3-PE-2.fa:2:30:7:1:true \\
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > 2_trim/{sample_name}.trim.log 2>&1")

cat(trim_cmd[1])

trim_dir <- "2_trim"
trim_dir %>% str_c("code/", .) %>% dir.create()

setCMD(trim_cmd, str_c("code/", trim_dir), 8, F)

# * * 2.3. QC for trimmed data --------------------------------------------

# reload trimmed data first

qc_dir <- "3_qc"
qc_dir %>% dir.create()

qc_cmd <- glue("fastqc -o {qc_dir} -t 8 0_fastq/{file_name} &")
cat(qc_cmd[1])

write.table(c("#!/bin/bash\n", qc_cmd), glue("code/{qc_dir}.sh"), quote = F, row.names = F, col.names = F)
# multiqc -o 3_qc -f -n qc.trim 3_qc/*.zip

# * * 2.4. Map ------------------------------------------------------------

ref <- "/lustre/user/liclab/lvyl/ref/hg19/refdata-cellranger-hg19-3.0.0/fasta/genome" # path to reference

map_dir <- "4_map"
map_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

map_cmd <- glue(
  "bowtie2 -p 20 --very-sensitive -X 2000 -x {ref} \\
  -1 {trim_dir}/{R1} -2 {trim_dir}/{R2} \\
  -S {map_dir}/{sample_name}.sam > {map_dir}/{sample_name}.log 2>&1")
cat(map_cmd[1])

setCMD(map_cmd, str_c("code/", map_dir), 6, T)

# * * 2.5. Sort -----------------------------------------------------------

sort_dir <- "5_sort"
sort_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

sort_cmd <- glue(
  "samtools view -@ 10 -bS {map_dir}/{sample_name}.sam | \\
  samtools sort -@ 10 > {sort_dir}/{sample_name}.bam")
cat(sort_cmd[1])

setCMD(sort_cmd, str_c("code/", sort_dir), 6, F)

# * * 2.6. Dedup ----------------------------------------------------------

picard <- "/lustre/user/liclab/lvyl/app/picard-2.18.21/picard.jar" # path to picard

dedup_dir <- "6_dedup"
dedup_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

dedup_cmd <- glue(
  "java -Xms2g -Xmx8g -XX:ParallelGCThreads=8 -jar {picard} MarkDuplicates \\
  I={sort_dir}/{sample_name}.bam \\
  O={dedup_dir}/{sample_name}.dedup.bam \\
  M={dedup_dir}/{sample_name}.dedup.txt \\
  REMOVE_DUPLICATES=true")
cat(dedup_cmd[1])

setCMD(dedup_cmd, str_c("code/", dedup_dir), 6, F)

# * * 2.7. Filter ---------------------------------------------------------

fil_dir <- "7_filter"
fil_dir %>% dir.create()

filter_cmd <- glue(
  "samtools view -h -f 2 -q 30 {dedup_dir}/{sample_name}.dedup.bam | \\
  egrep -v 'MT|GL|JH' | samtools sort -@ 4 -O bam > {fil_dir}/{sample_name}.filter.bam &")
cat(filter_cmd[1])

write.table(c("#!/bin/bash\n", filter_cmd), glue("code/{fil_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.8. Finger ---------------------------------------------------------

fin_dir <- "8_finger"
fin_dir %>% dir.create()

index_cmd <- glue("samtools index {fil_dir}/{sample_name}.filter.bam {fil_dir}/{sample_name}.filter.bam.bai &")
cat(index_cmd)

write.table(c("#!/bin/bash\n", index_cmd), "code/8_index.sh", quote = F, row.names = F, col.names = F)

finger_group <- tapply(sample_name, gsub(".*-", "", sample_name), function(x) glue("{fil_dir}/{x}.filter.bam"))
finger_group %<>% map(~ paste(.x, collapse = " "))

finger_cmd <- glue("
  plotFingerprint -b {finger_group} -o {fin_dir}/{names(finger_group)}.finger.pdf &")
cat(finger_cmd[1])

write.table(c("#!/bin/bash\n", finger_cmd), glue("code/{fin_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.9. Peak -----------------------------------------------------------

peak_dir <- "9_peak"
peak_dir %>% dir.create()

input_sample <- grep("in", sample_name, value = T)
chip_sample <- grep("in", sample_name, value = T, invert = T)
input_sample %<>% rep(each = 3)

org <- "hs"

peak_cmd <- glue(
  "macs2 callpeak -t {fil_dir}/{chip_sample}.filter.bam -c {fil_dir}/{input_sample}.filter.bam \\
  -f BAMPE -g {org} --broad --keep-dup all -n {chip_sample} \\
  --outdir {peak_dir} > {peak_dir}/{chip_sample}.log 2>&1 &")
cat(peak_cmd[1])

write.table(c("#!/bin/bash\n", peak_cmd), glue("code/{peak_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.10. Bigwig --------------------------------------------------------

bw_dir <- "10_bw"
bw_dir %>% dir.create()

bw_cmd <- glue("
  bamCoverage --normalizeUsing RPKM -of bigwig \\
  -b {fil_dir}/{sample_name}.filter.bam \\
  -o {bw_dir}/{sample_name}.bw > {bw_dir}/{sample_name}.log 2>&1 &")
cat(bw_cmd[1])

bw_comp_cmd <- glue(
  "bamCompare -p 2 -of bigwig --operation ratio \\
  -b1 {fil_dir}/{chip_sample}.filter.bam \\
  -b2 {fil_dir}/{input_sample}.filter.bam \\
  -o {bw_dir}/{chip_sample}.comp.bw > {bw_dir}/{chip_sample}.comp.log 2>&1 &")
cat(bw_comp_cmd[1])

write.table(c("#!/bin/bash\n", bw_cmd), glue("code/{bw_dir}.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", bw_comp_cmd), glue("code/{bw_dir}_comp.sh"), quote = F, row.names = F, col.names = F)

# * * 2.11. Heat ----------------------------------------------------------

chip_type <- chip_sample %>% str_replace(".*-", "")
bw_group <- tapply(chip_sample, chip_type, c) %>%
  map(~ glue("{bw_dir}/{.x}.comp.bw")) %>%
  map(~ str_c(.x, collapse = " "))

heat_dir <- "11_heat"
heat_dir %>% dir.create()

bed_file <- "all.bed"
out_prefix <- str_c("all_", names(bw_group))

mtx_cmd <- glue(
  "computeMatrix scale-regions -p 6 -R {heat_dir}/{bed_file} \\
  -S {bw_group} -o {heat_dir}/{out_prefix}_mtx.gz &"
)
cat(mtx_cmd[1])

heat_cmd <- glue(
  "plotHeatmap -m {heat_dir}/{out_prefix}_mtx.gz --colorMap RdBu_r \\
  -o {heat_dir}/{out_prefix}_heatmap.pdf \\
  --outFileNameMatrix {heat_dir}/{out_prefix}_value.txt \\
  --outFileSortedRegions {heat_dir}/{out_prefix}_region.bed &"
)
cat(heat_cmd[1])

write.table(c("#!/bin/bash\n", mtx_cmd), glue("code/{heat_dir}_1.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", heat_cmd), glue("code/{heat_dir}_2.sh"), quote = F, row.names = F, col.names = F)

# * 3. Function -----------------------------------------------------------

setCMD <- function(cmd, dir = "", sepN = 1, clu = F) {
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

