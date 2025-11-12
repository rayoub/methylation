
library(randomForest)
library(glmnet)

source("common.R")
source("mcf.R")

rf_model <- loadSavedModel(REF_GSE_ID)

# **** calibration model ****

anno <- loadSavedAnno(REF_GSE_ID)
betas <- loadSavedBetas(REF_GSE_ID)

scores <- predict(rf_model, betas, type="prob")
rownames(scores) <- sub("_.*", "", rownames(scores))
y <- anno$`methylation class:ch1`

suppressWarnings(
	calfit <- cv.glmnet(
					y=y,
					x=scores,
					family="multinomial",
					type.measure="mse",
					alpha=0,
					nlambda=100,
					lambda.min.ratio=10^-6,
					parallel=TRUE
				)
)

saveRDS(calfit, file=file.path("results", paste0(REF_GSE_ID,"_calfit")))

# **** prediction (with score calibration) ****

betas <- loadSavedBetas(VAL_GSE_ID)
scores <- predict(rf_model, betas, type="prob")
rownames(scores) <- sub("_.*","",rownames(scores))
scores <- predict(calfit$glmnet.fit,
						newx=scores,
						type="response",
						s=calfit$lambda.1se)[,,1] # use lambda estimated by 10fold CVlambda

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



