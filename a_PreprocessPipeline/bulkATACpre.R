
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
  fastqcDir = "fastqc",
  trimmomaticDir = "app/Trimmomatic-0.39/trimmomatic-0.39.jar",
  trimAdapter = "app/Trimmomatic-0.39/adapters/NexteraPE-PE.fa",
  bt2Dir = "bowtie2",
  mapRef = "fasta/genome",
  stDir = "samtools",
  picardDir = "app/picard.jar",
  bedtDir = "bedtools",
  macs2Dir = "macs2",
  dtDir = "",
  bdg2bwDir = "app/bedGraphToBigWig",
  refLen = "fasta/genome.fa.fai",
  idrDir = "idr",
  homerDir = "",
  knownMotif = "app/homer/data/knownTFs/vertebrates/known.motifs"
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

qc <- settingList$fastqcDir

qc_cmd <- glue("{qc} -o {qc_dir} -t 8 0_fastq/{file_name} &")
cat(qc_cmd[1])

write.table(c("#!/bin/bash\n", qc_cmd), glue("code/{qc_dir}.sh"), quote = F, row.names = F, col.names = F)
# multiqc -o 1_qc -f -n qc.raw 1_qc/*.zip

# * * 2.2. Trim -----------------------------------------------------------

dir.create("2_trim")
dir.create(".2_untrim")

trim <- settingList$trimmomaticDir
adapter <- settingList$trimAdapter

trim_cmd <- glue(
  "java -jar {trim} PE -threads 10 \\
  0_fastq/{R1} 0_fastq/{R2} \\
  2_trim/{sample_name}_1.trim.fastq.gz \\
  .2_untrim/{sample_name}_1.untrim.fastq.gz \\
  2_trim/{sample_name}_2.trim.fastq.gz \\
  .2_untrim/{sample_name}_2.untrim.fastq.gz \\
  ILLUMINACLIP:{adapter}:2:30:7:1:true \\
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

bt2 <- settingList$bt2Dir
ref <- settingList$mapRef

map_dir <- "4_map"
map_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

map_cmd <- glue(
  "{bt2} -p 20 --very-sensitive -X 2000 -x {ref} \\
  -1 {trim_dir}/{R1} -2 {trim_dir}/{R2} \\
  -S {map_dir}/{sample_name}.sam > {map_dir}/{sample_name}.log 2>&1")
cat(map_cmd[1])

setCMD(map_cmd, str_c("code/", map_dir), 1, T)
# multiqc -o 4_map -f -n map 4_map/*.log

# * * 2.5. Sort -----------------------------------------------------------

st <- settingList$stDir

sort_dir <- "5_sort"
sort_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

sort_cmd <- glue(
  "{st} view -@ 10 -bS {map_dir}/{sample_name}.sam | \\
  {st} sort -@ 10 > {sort_dir}/{sample_name}.bam")
cat(sort_cmd[1])

setCMD(sort_cmd, str_c("code/", sort_dir), 1, T)

# * * 2.6. Dedup ----------------------------------------------------------

picard <- settingList$picardDir

dedup_dir <- "6_dedup"
dedup_dir %T>% dir.create() %>% str_c("code/", .) %>% dir.create()

dedup_cmd <- glue(
  "java -Xms2g -Xmx8g -XX:ParallelGCThreads=8 -jar {picard} MarkDuplicates \\
  I={sort_dir}/{sample_name}.bam \\
  O={dedup_dir}/{sample_name}.dedup.bam \\
  M={dedup_dir}/{sample_name}.dedup.txt \\
  REMOVE_DUPLICATES=true > {dedup_dir}/{sample_name}.dedup.log 2>&1")
cat(dedup_cmd[1])

setCMD(dedup_cmd, str_c("code/", dedup_dir), 1, T)

# * * 2.7. Filter ---------------------------------------------------------

fil_dir <- "7_filter"
fil_dir %>% dir.create()

filter_cmd <- glue(
  "{st} view -h -f 2 -q 30 {dedup_dir}/{sample_name}.dedup.bam | \\
  egrep -v '\\bMT|\\bGL|\\bJH' | {st} sort -@ 4 -O bam > {fil_dir}/{sample_name}.filter.bam &")
cat(filter_cmd[1])

write.table(c("#!/bin/bash\n", filter_cmd), glue("code/{fil_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.8. Insert ---------------------------------------------------------

picard <- settingList$picardDir

ins_dir <- "8_insert"
ins_dir %>% dir.create()

ins_cmd <- glue(
  "java -Xms2g -Xmx8g -XX:ParallelGCThreads=8 -jar {picard} CollectInsertSizeMetrics \\
  I={fil_dir}/{sample_name}.filter.bam \\
  O={ins_dir}/{sample_name}.txt \\
  H={ins_dir}/{sample_name}.pdf > {ins_dir}/{sample_name}.log 2>&1 &")
cat(ins_cmd[1])

write.table(c("#!/bin/bash\n", ins_cmd), glue("code/{ins_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.9. Shift ----------------------------------------------------------

bedt <- settingList$bedtDir

shift_dir <- "9_shift"
shift_dir %>% dir.create()

shift_cmd <- glue(
  "{bedt} bamtobed -i {fil_dir}/{sample_name}.filter.bam | {awk_cmd} \\
  > {shift_dir}/{sample_name}.shift.bed &",
  awk_cmd = "awk -F '\\t' 'BEGIN {OFS = FS}{ if ($6 == \"+\") {$2 = $2 + 4} else if ($6 == \"-\") {$3 = $3 - 5} print $0}'")
cat(shift_cmd[1])

write.table(c("#!/bin/bash\n", shift_cmd), glue("code/{shift_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.10. Peak ----------------------------------------------------------

macs2 <- settingList$macs2Dir

peak_dir <- "10_peak"
peak_dir %>% dir.create()

org <- "hs"
org <- "mm"

peak_cmd <- glue("
  {macs2} callpeak -t {shift_dir}/{sample_name}.shift.bed \\
  -g {org} --nomodel --shift -100 --extsize 200 --keep-dup all -n {sample_name} \\
  -B --outdir {peak_dir} > {peak_dir}/{sample_name}.log 2>&1 &")
cat(peak_cmd[1])

group_name <- tapply(sample_name, str_sub(sample_name, 1, -3), c)

group_peak_cmd <- glue(
  "{macs2} callpeak -t {peak_group} \\
  -g {org} --nomodel --shift -100 --extsize 200 --keep-dup all --SPMR -n group_{group} \\
  -B --outdir {peak_dir} > {peak_dir}/group_{group}.log 2>&1 &",
  peak_group = map_chr(group_name, ~ paste(glue("{shift_dir}/{.x}.shift.bed"), collapse = " ")),
  group = names(group_name))
cat(group_peak_cmd[1])

write.table(c("#!/bin/bash\n", peak_cmd), glue("code/{peak_dir}.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", group_peak_cmd), glue("code/{peak_dir}_group.sh"), quote = F, row.names = F, col.names = F)

# * * 2.11. Bigwig --------------------------------------------------------

bw_dir <- "11_bw"
bw_dir %>% dir.create()

dt <- settingList$dtDir
bdg2bw <- settingList$bdg2bwDir
ref_len <- settingList$refLen

index_cmd <- glue("{st} index {fil_dir}/{sample_name}.filter.bam &")
cat(index_cmd[1])

bw_cmd <- glue("
  {dt}bamCoverage --normalizeUsing RPKM -of bigwig -b {fil_dir}/{sample_name}.filter.bam \\
  -o {bw_dir}/{sample_name}.bw &")
cat(bw_cmd[1])

group_bw_cmd <- glue(
  "{bdg2bw} {peak_dir}/group_{group}_treat_pileup.bdg {ref_len} {bw_dir}/{group}.bw &",
  group = names(group_name))
cat(group_bw_cmd[1])

sample_name <- list.files(bw_dir, "-[0-9].bw") %>% str_replace(".bw$", "")
group_name <- tapply(sample_name, str_sub(sample_name, 1, -3), c)

cor_cmd <- glue(
  "{dt}multiBigwigSummary bins -p 20 -b {bws} -o {bw_dir}/bw_cor.npz &",
  bws = str_c(glue("{bw_dir}/{sample_name}.bw"), collapse = " \\\n"))
cat(cor_cmd)

group_cor_cmd <- glue(
  "{dt}multiBigwigSummary bins -p 20 -b {group_bws} -o {bw_dir}/group_bw_cor.npz &",
  group_bws = str_c(glue("{bw_dir}/{names(group_name)}.bw"), collapse = " \\\n"))
cat(group_cor_cmd)

corplot_cmd <- glue(
  "{dt}plotCorrelation -in {bw_dir}/bw_cor.npz -c pearson -p heatmap --plotNumbers -o {bw_dir}/cor_heatmap.pdf --colorMap RdBu_r &")
cat(corplot_cmd)

group_corplot_cmd <- glue(
  "{dt}plotCorrelation -in {bw_dir}/group_bw_cor.npz -c pearson -p heatmap --plotNumbers -o {bw_dir}/group_cor_heatmap.pdf --colorMap RdBu_r &")
cat(group_corplot_cmd)

write.table(c("#!/bin/bash\n", index_cmd), glue("code/{bw_dir}_index.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", bw_cmd), glue("code/{bw_dir}.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", group_bw_cmd), glue("code/{bw_dir}_group.sh"), quote = F, row.names = F, col.names = F)

write.table(c("#!/bin/bash\n", cor_cmd), glue("code/{bw_dir}_cor.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", group_cor_cmd), glue("code/{bw_dir}_cor_group.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", corplot_cmd), glue("code/{bw_dir}_plot_cor.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", group_corplot_cmd), glue("code/{bw_dir}_plot_cor_group.sh"), quote = F, row.names = F, col.names = F)

# * * 2.12. IDR -----------------------------------------------------------

idr_dir <- "12_idr"
idr_dir %>% dir.create()

idr <- settingList$idrDir

idr_cmd <- map(group_name, ~ as_tibble(combn(.x, 2)) %>% map(~ paste(glue("{peak_dir}/{.x}_peaks.narrowPeak"), collapse = " "))) %>% 
  imap(~ glue("{idr} --samples {.x} --peak-list {peak_dir}/group_{.y}_peaks.narrowPeak \\
              --input-file-type narrowPeak --output-file {idr_dir}/{.y}")) %>% 
  map(~ imap(.x, ~ glue("{.x}-{.y}.idr.narrowPeak &"))) %>% unlist() %>% unname()
cat(idr_cmd[1])

write.table(c("#!/bin/bash\n", idr_cmd), glue("code/{idr_dir}.sh"), quote = F, row.names = F, col.names = F)

group_index <- map(group_name, ~ choose(length(.x), 2)) %>% imap(~ rep(.y, .x)) %>% unlist() %>% unname()
idr_file <- list.files(idr_dir, "idr", full.names = T)
idr_data <- map(idr_file, ~ read_delim(.x, "\t", col_names = F))
idr_peak <- map(idr_data, ~ .x[.x[[5]] > 540, ]) %>%
  map(~ str_c(.x[[1]], .x[[2]], .x[[3]], sep = "_")) %>%
  tapply(group_index, c) %>%
  map(~ unlist(.x) %>% unique)

group_file <- list.files(peak_dir, "^group.*Peak", full.names = T)
group_data <- map(group_file, ~ read.table(.x, sep = "\t", stringsAsFactors = F))
group_peak <- map(group_data, ~ str_c(.x[[1]], .x[[2]], .x[[3]], sep = "_"))

groupData <- pmap(list(idr_peak, group_peak, group_data), function(a, b, c) {c[b %in% a, ]})

iwalk(group_data, ~ write.table(.x, glue("{idr_dir}/{.y}.narrowPeak"), col.names = F, row.names = F, quote = F, sep = "\t"))

# * * 2.13. Motif ---------------------------------------------------------

motif_dir <- "13_motif"
motif_dir %>% dir.create()
region_file <- list.files(motif_dir, "bed")

homer <- settingList$homerDir

org <- "hg19"
org <- "mm10"

motif_cmd <- glue(
  "{homer}findMotifsGenome.pl {motif_dir}/{region_file} {org} {motif_dir}/{rigion_name} \\
  -size 200 -len 8,10,12 > {region_file}/{region_name}.log 2>&1 &",
  region_name = str_replace(region_file, ".bed", ""))
cat(motif_cmd[1])

write.table(c("#!/bin/bash\n", motif_cmd), glue("code/{motif_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.14. Heat ----------------------------------------------------------

heat_dir <- "14_heat"
heat_dir %>% dir.create()

sample_name <- list.files(bw_dir, "-[123].bw") %>% str_replace(".bw$", "")
group_name <- tapply(sample_name, str_sub(sample_name, 1, -3), c)

mat_cmd <- glue(
  "{dt}computeMatrix reference-point -a 2000 -b 2000 -p 20 -R {heat_dir}/hg19geneTSS.bed -S \\
  {bws} -o {heat_dir}/TSS_mtx.gz &",
  bws = str_c(glue("{bw_dir}/{names(group_name)}.bw"), collapse = " \\\n"))
cat(mat_cmd)

heat_cmd <- glue(
  "{dt}plotHeatmap -m {heat_dir}/TSS_mtx.gz --colorMap RdBu_r -o {heat_dir}/TSS.heatmap.pdf \\
  --outFileNameMatrix {heat_dir}/TSS_value.txt --outFileSortedRegions {heat_dir}/TSS_region.bed &")
cat(heat_cmd)

write.table(c("#!/bin/bash\n", mat_cmd), glue("code/{heat_dir}_mtx.sh"), quote = F, row.names = F, col.names = F)
write.table(c("#!/bin/bash\n", heat_cmd), glue("code/{heat_dir}.sh"), quote = F, row.names = F, col.names = F)

# * * 2.15. Anno ----------------------------------------------------------

anno_dir <- "15_anno"
anno_dir %>% dir.create()

known_motif <- settingList$knownMotif

peak_file <- list.files(anno_dir, ".homer")

anno_cmd <- glue(
  "{homer}annotatePeaks.pl {anno_dir}/{peak_file} {org} -m {known_motif} -mscore > {anno_dir}/{peak_name}.txt &",
  peak_name = str_remove(peak_file, ".homer")
)
cat(anno_cmd)

write.table(c("#!/bin/bash\n", anno_cmd), glue("code/{anno_dir}.sh"), quote = F, row.names = F, col.names = F)

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
