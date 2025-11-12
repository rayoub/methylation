
library(randomForest)

source("common.R")
source("mcf.R")

# **** prediction ****

rf_model <- loadSavedModel(REF_GSE_ID)
betas <- loadSavedBetas(VAL_GSE_ID)
scores <- predict(rf_model, betas, type="prob")
rownames(p) <- sub("_.*","",rownames(scores))

# **** evaluation ****

anno <- loadSavedAnno(VAL_GSE_ID)

anno <- anno[,"methylation class:ch1", drop=FALSE]
colnames(anno) <- c("MC")

mcf <- unlist(ifelse(anno[, "MC"] %in% names(MCF_LOOKUP), MCF_LOOKUP[anno[, "MC"]], anno[, "MC"]))
mc_pred <- colnames(scores)[max.col(scores)]
anno <- cbind(anno, MCF = mcf, MC_PRED = mc_pred)

mcf_pred <- unlist(ifelse(anno[, "MC_PRED"] %in% names(MCF_LOOKUP), MCF_LOOKUP[anno[, "MC_PRED"]], anno[, "MC_PRED"]))
anno <- cbind(anno, MCF_PRED = mcf_pred)

sum(anno$MC == anno$MC_PRED) / nrow(anno)
sum(anno$MCF == anno$MCF_PRED) / nrow(anno)

# descr								mc_pred			mcf_pred
# without batch adjustments			0.8577899		0.928442
# with batch adjustments			0.8804348		0.9347826
