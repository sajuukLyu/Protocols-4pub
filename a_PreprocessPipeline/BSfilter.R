
# MESSAGE -----------------------------------------------------------------
# 
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
# 
# require: data.table
# 
# usage: Rscript BSfilter.R -i [input] -d [DSS output] -b [BDG output] -p [non-convert rate]
# 
# ---

# load package ------------------------------------------------------------

options(scipen = 200)

suppressMessages(library(data.table))

message("LOG\t", format(Sys.time(), "%H:%M:%S"), "\tpackage loaded")

# solve args --------------------------------------------------------------

args <- commandArgs(T)
INPUT <- args[which(args == "-i") + 1]
DSS <- args[which(args == "-d") + 1]
BDG <- args[which(args == "-b") + 1]
P <- as.numeric(args[which(args == "-p") + 1])

message("ARG\tinput:\t", INPUT,
        "\n\tDSSout:\t", DSS,
        "\n\tBDGout:\t", BDG,
        "\n\tp:\t", P)

# analyze -----------------------------------------------------------------

message("LOG\t", format(Sys.time(), "%H:%M:%S"), "\tloading data ...")
data <- fread(INPUT)
colnames(data) <- c("chr", "loci", "chain", "m", "u", "type", "context")

message("LOG\t", format(Sys.time(), "%H:%M:%S"), "\tfiltering ...")
data <- data[!grepl("^(GL|MT)", chr)]

message("LOG\t", format(Sys.time(), "%H:%M:%S"), "\ttesting ...")
data[, id := paste(m, u, sep = "_")]
setkey(data, id)

meth_test <- function(m, a, p) {
  ifelse(binom.test(m, a, p, "greater")$p.value >= 0.01, "u", "m")
}

data[, res := meth_test(.SD[1][1, m], sum(.SD[1]), P), by = id, .SDcols = c("m", "u")]
data[res == "u", u := m + u][res == "u", m := 0]

message("LOG\t", format(Sys.time(), "%H:%M:%S"), "\tsorting ...")
setkey(data, chr, loci)

# output data -------------------------------------------------------------

data[, value := m / (m + u)]

message("LOG\t", format(Sys.time(), "%H:%M:%S"), "\toutputing ...")
fwrite(data[, .(chr, loci, m + u, m)], DSS, sep = "\t", col.names = F)
fwrite(data[, .(chr, loci, loci + 1, value)], BDG, sep = "\t", col.names = F)

message("LOG\t", format(Sys.time(), "%H:%M:%S"), "\done!")
