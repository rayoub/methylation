
library(ranger)
library(openxlsx)
library(glmnet)
library(here)

source(here("R","globals.R"))
source(here("R","loading.R"))
source(here("R","prediction.R"))

id <- "NIH_EPIC"
betas <- loadSavedBetas(id)

scores <- predictSampleScores(betas)

res <- evaluateDiagnosticSampleScores(scores)
#res <- evaluateGEOSampleScores(id, scores)

sum(res$MC == res$MC_PRED) / nrow(res)
sum(res$MCF == res$MCF_PRED) / nrow(res)

write.xlsx(res, file="nih_epic.xlsx", rowNames=TRUE)

# old randomForest
# descr															mc_pred			mcf_pred
# without batch adjustments and 30000/10000 filter				0.8577899		0.928442
# + batch adjustments											0.8804348		0.9347826
# + calibration													0.8405797 		0.9211957*

# switching to ranger
# descr															mc_pred			mcf_pred
# 100000/10000 filter + batch adjustments						0.8894928		0.9447464
# no filters/10000 + batch adjustments							0.8958333		0.9456522
# no filters/30000 + batch adjustments							0.8967391		0.9483696

#																0.8967391		0.9474638
#																0.9076087		0.9547101



