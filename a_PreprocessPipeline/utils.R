
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
