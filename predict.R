
library(randomForest)

source("common.R")
source("mcf.R")

rf.model <- loadSavedModel(REF_GSE_ID)
betas <- loadSavedBetas(VAL_GSE_ID)
p <- predict(rf.model, betas, type="prob")
rownames(p) <- sub("_.*","",rownames(p))

anno <- loadSavedAnno(VAL_GSE_ID)

anno <- anno[,"methylation class:ch1", drop=FALSE]
colnames(anno) <- c("MC")

mcf <- unlist(ifelse(anno[, "MC"] %in% names(MCF_LOOKUP), MCF_LOOKUP[anno[, "MC"]], anno[, "MC"]))
mc_pred <- colnames(p)[max.col(p)]
anno <- cbind(anno, MCF = mcf, MC_PRED = mc_pred)

mcf_pred <- unlist(ifelse(anno[, "MC_PRED"] %in% names(MCF_LOOKUP), MCF_LOOKUP[anno[, "MC_PRED"]], anno[, "MC_PRED"]))
anno <- cbind(anno, MCF_PRED = mcf_pred)
