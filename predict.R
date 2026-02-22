
library(dplyr)
library(tidyr)
library(tibble)
library(openxlsx)
library(here)

source(here("R","constants.R"))
source(here("R","prediction.R"))

# for lab samples
batch_id <- "BATCH1"
scores <- predictLabSampleScores(batch_id)
res <- evaluateLabSampleScores(scores)

# # for geo samples
# gse_id <- REF_GSE_ID
# scores <- predictGeoSampleScores(batch_id)
# res <- evaluateGeoSampleScores(scores)
# 
# # geo evaluation
# sum(res$MC == res$MC_PRED) / nrow(res)
# sum(res$MCF == res$MCF_PRED) / nrow(res)
# 
# # write geo results
# write.xlsx(res, file=here(paste0(gse_id, ".xlsx")), rowNames=TRUE)
