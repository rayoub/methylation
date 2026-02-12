
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)
library(here)

source(here("R","constants.R"))
source(here("R","loading.R"))
source(here("R","prediction.R"))


id <- "DIAG3"
scores <- predictSampleScores(id)

res <- evaluateDiagnosticSampleScores(scores)
#res <- evaluateGEOSampleScores(id, scores)

#sum(res$MC == res$MC_PRED) / nrow(res)
#sum(res$MCF == res$MCF_PRED) / nrow(res)

write.xlsx(res, file=here(paste0(id, ".xlsx")), rowNames=TRUE)
