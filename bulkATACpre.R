
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

trim <- "/lustre/user/liclab/lvyl/app/Trimmomatic-0.39"
adapter <- "NexteraPE-PE.fa"

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

setCMD(trim_cmd, str_c("code/", trim_dir), 7, F)

# * * 2.3. QC for trimmed data --------------------------------------------

# reload trimmed data first

qc_dir <- "3_qc"
qc_dir %>% dir.create()

qc_cmd <- glue("fastqc -o {qc_dir} -t 8 {trim_dir}/{file_name} &")
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

# * * 2.8. Insert ---------------------------------------------------------

picard <- "/lustre/user/liclab/lvyl/app/picard-2.18.21/picard.jar" # path to picard

ins_dir <- "8_insert"
ins_dir %>% dir.create()

ins_cmd <- glue(
  "java -Xms2g -Xmx8g -XX:ParallelGCThreads=8 -jar {picard} CollectInsertSizeMetrics \\
  I={fil_dir}/{sample_name}.filter.bam \\
  O={ins_dir}/{sample_name}.txt \\
  H={ins_dir}/{sample_name}.pdf > {ins_dir}/{sample_name}.log 2>&1")
cat(ins_cmd[1])

write.table(c("#!/bin/bash\n", ins_cmd), glue("code/{ins_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.9. Shift ----------------------------------------------------------

shift_dir <- "9_shift"
shift_dir %>% dir.create()

shift_cmd <- glue(
  "bedtools bamtobed -i {fil_dir}/{sample_name}.filter.bam | {awk_cmd} \\
  > {shift_dir}/{sample_name}.shift.bed &",
  awk_cmd = "awk -F '\\t' 'BEGIN {OFS = FS}{ if ($6 == \"+\") {$2 = $2 + 4} else if ($6 == \"-\") {$3 = $3 - 5} print $0}'")
cat(shift_cmd[1])

write.table(c("#!/bin/bash\n", shift_cmd), glue("code/{shift_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.10. Peak ----------------------------------------------------------

peak_dir <- "10_peak"
peak_dir %>% dir.create()

org <- "hs"
org <- "mm"

peak_cmd <- glue("
  macs2 callpeak -t {shift_dir}/{sample_name}.shift.bed \\
  -g {org} --nomodel --shift -100 --extsize 200 --keep-dup all -n {sample_name} \\
  -B --outdir {peak_dir} > {peak_dir}/{sample_name}.log 2>&1 &")
cat(peak_cmd[1])

group_name <- list(
  c("ENR", "XL5"),
  c("NC", "XL4"),
  c("VE", "XL1")
)
names(group_name) <- c("ENR", "NC", "VE")

peak_group <- map_chr(group_name, ~ paste(glue("{shift_dir}/{.x}.shift.bed"), collapse = " "))

group_peak_cmd <- glue(
  "macs2 callpeak -t {peak_group} \\
  -g {org} --nomodel --shift -100 --extsize 200 --keep-dup all --SPMR -n group_{group} \\
  -B --outdir {peak_dir} > {peak_dir}/group_{group}.log 2>&1 &",
  group = names(group_name))
cat(group_peak_cmd[1])

write.table(c("#!/bin/bash\n", peak_cmd), glue("code/{peak_dir}.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", group_peak_cmd), glue("code/{peak_dir}_group.sh"), quote = F, row.names = F, col.names = F)

# * * 2.11. Bigwig --------------------------------------------------------

bw_dir <- "11_bw"
bw_dir %>% dir.create()

index_cmd <- glue("samtools index {fil_dir}/{sample_name}.filter.bam &")
cat(index_cmd[1])

bw_cmd <- glue("
  bamCoverage --normalizeUsing RPKM -of bigwig -b {fil_dir}/{sample_name}.filter.bam \\
  -o {bw_dir}/{sample_name}.bw &")
cat(bw_cmd[1])

write.table(c("#!/bin/bash\n", index_cmd), glue("code/{bw_dir}_index.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", bw_cmd), glue("code/{bw_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.12. IDR -----------------------------------------------------------

map(group_name, ~ combn(.x, 2))

idr_cmd <- map(group_name, ~ as.data.frame(combn(.x, 2)) %>% map(~ paste(glue("10_peak/{.x}_peaks.narrowPeak"), collapse = " "))) %>% 
  imap(~ glue("idr --samples {.x} --peak-list 10_peak/group_{.y}_peaks.narrowPeak \\
              --input-file-type narrowPeak --output-file 12_idr/{.y}")) %>% 
  map(~ imap(.x, ~ glue("{.x}-{.y}.idr.narrowPeak &"))) %>% unlist() %>% unname()
cat(idr_cmd[1])

dir.create("12_idr")

write.table(c("#!/bin/bash\n", idr_cmd), "code/12_idr.sh", quote = F, row.names = F, col.names = F)

idrFile <- list.files("12_idr", full.names = T)
idrPeak <- map(idrFile, ~ read.table(.x, sep = "\t", stringsAsFactors = F, comment.char = ""))
names(idrPeak) <- names(group_name)

map_int(idrPeak, nrow)

idrPeak %<>% map(~ .x[.x[[5]] > 540, ])

iwalk(idrPeak, ~ write.table(.x, glue("12_idr/{.y}.narrowPeak"), col.names = F, row.names = F, quote = F, sep = "\t"))

# * * 2.13. Motif ---------------------------------------------------------

dir.create("13_motif")
region_file <- list.files("13_motif", "bed")

motif_cmd <- paste0("findMotifsGenome.pl 13_motif/", region_file, " mm10 13_motif/",
                    gsub(".bed", "", region_file), " -size 200 -len 8,10,12 > 13_motif/",
                    gsub(".bed", "", region_file), ".log 2>&1 &")
cat(motif_cmd[1])

write.table(c("#!/bin/bash\n", motif_cmd), "code/13_motif.sh", quote = F, row.names = F, col.names = F)

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
