
# MESSAGE -----------------------------------------------------------------
#
# author: Yulin Lyu
# email: lvyulin@pku.edu.cn
#
# require: R whatever
#
# ---

# * 1. Load packages ------------------------------------------------------

setwd("exampleData/ATAC")

# grammar
library(tidyverse)
library(magrittr)
library(glue)
library(data.table)

# analysis
library(ggplot2)
library(ggsci)
library(scales)
library(patchwork)

# * 2. TSS enrichment -----------------------------------------------------

dir.create("graphics")

usedSample <- c("F-H1", "P-H1", "XF-H1")

TSSdata <- fread("qc/TSS_value.txt", skip = 3) %>% as.matrix()
TSSdata[is.nan(TSSdata)] <- 0

TSSmean <- apply(TSSdata, 2, mean)
TSSsam <- split(TSSmean, rep(usedSample, each = 400))

enrichScore <- map_dbl(TSSsam, ~ {sum(.x[191:210]) / sum(.x[c(1:10, 391:400)])})
saveRDS(enrichScore, "enrichScore.rds")

plotData <- data.table(
  pos = rep(1:400, 3),
  val = TSSmean,
  sample = rep(usedSample, each = 400)
)

usedLabel <- str_c(usedSample, ": ", round(enrichScore, 2))
usedColor <- ArchR::ArchRPalettes$stallion[seq_along(usedSample)]
names(usedColor) <- usedSample

ggplot(plotData, aes(x = pos, y = val)) +
  geom_line(aes(color = sample)) +
  scale_color_manual(values = usedColor, labels = usedLabel) +
  scale_x_continuous(
    breaks = c(0, 200, 400),
    labels = c("-2k", "TSS", "+2k")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Mean signal") +
  theme(
    aspect.ratio = 0.9,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0, 1),
    legend.justification = c("left", "top")
  )

ggsave("graphics/TSSenrich.png", width = 4, height = 4)

# * 3. Fragment length ----------------------------------------------------

fragFile <- list.files("qc", "-[123].txt", full.names = T)
fragSample <- str_replace_all(fragFile, ".*/|.txt", "")

fragList <- map(fragFile, ~ fread(.x, skip = 10)) %>% set_names(fragSample)

fragData <- fragList %>% imap(~ {.x[, sample := .y]}) %>% reduce(rbind)
fragData[, type := str_remove(sample, "-.$")][]
colnames(fragData)[1:2] <- c("pos", "count")

plotData <- map(usedSample, ~ {fragData[type == .x]})

usedColor <- ArchR::ArchRPalettes$stallion %>% unname()

plotList <- map(plotData, ~ {
  ggplot(.x, aes(x = pos, y = count)) +
    geom_line(aes(color = sample)) +
    scale_color_manual(values = usedColor) +
    labs(x = "", y = "") +
    theme(
      aspect.ratio = 1,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.position = c(1, 1),
      legend.justification = c("right", "top")
    )
})

reduce(plotList, `+`) + plot_layout(ncol = 1)

ggsave("graphics/fragLength.png", width = 4, height = 3 * length(usedSample))

