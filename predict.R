
library(ranger)
library(openxlsx)
library(glmnet)


source("common.R")
source("mcf.R")
source("mc.R")

do.predict <- function(rf_model, betas) {

	rf_vars <- rf_model$forest$independent.variable.names
	filtered_betas <- betas[match(rf_vars, colnames(betas))]
	p <- predict(rf_model, data = filtered_betas)
	scores <- p$predictions
	rownames(scores) <- sub("_.*", "", rownames(filtered_betas))

	return(scores)
}

do.evaluateGEO <- function (gse_id, scores) {

	anno <- loadSavedAnno(gse_id)
	if ("methylation class:ch1" %in% names(anno)) {

		res <- data.frame(MC = anno[,"methylation class:ch1"])

		mc_pred <- colnames(scores)[max.col(scores)]
		mcf <- unlist(ifelse(res[, "MC"] %in% names(MCF_LOOKUP), MCF_LOOKUP[res[, "MC"]], res[, "MC"]))
		res <- cbind(res, MC_PRED = mc_pred, MCF = mcf)

		mcf_pred <- unlist(ifelse(res[, "MC_PRED"] %in% names(MCF_LOOKUP), MCF_LOOKUP[res[, "MC_PRED"]], res[, "MC_PRED"]))
		res <- cbind(res, MCF_PRED = mcf_pred)

		mc_descr <- ifelse(res[, "MC"] %in% names(MC_DESCR), MC_DESCR[res[, "MC"]], NA)
		res <- cbind(res, MC_DESCR = mc_descr)
	}
	else {

		samples <- rownames(anno)	
		mc_pred <- colnames(scores)[max.col(scores)]
		res <- data.frame(SAMPLE_ID = samples, MC_PRED = mc_pred)

		mcf_pred <- unlist(ifelse(res[, "MC_PRED"] %in% names(MCF_LOOKUP), MCF_LOOKUP[res[, "MC_PRED"]], res[, "MC_PRED"]))
		res <- cbind(res, MCF_PRED = mcf_pred)

		mc_descr <- ifelse(res[, "MC_PRED"] %in% names(MC_DESCR), MC_DESCR[res[, "MC_PRED"]], NA)
		res <- cbind(res, MC_DESCR = mc_descr)

		nm_descr <- anno$`tumor type/grade:ch1`
		res <- cbind(res, NM_DESCR = nm_descr)
	}

	return (res)
}

do.evaluateDIAG <- function (scores) {

	samples <- rownames(scores)	
	mc_pred <- colnames(scores)[max.col(scores)]
	res <- data.frame(SAMPLE_ID = samples, MC_PRED = mc_pred)

	mcf_pred <- unlist(ifelse(res[, "MC_PRED"] %in% names(MCF_LOOKUP), MCF_LOOKUP[res[, "MC_PRED"]], res[, "MC_PRED"]))
	res <- cbind(res, MCF_PRED = mcf_pred)

	mc_descr <- ifelse(res[, "MC_PRED"] %in% names(MC_DESCR), MC_DESCR[res[, "MC_PRED"]], NA)
	res <- cbind(res, MC_DESCR = mc_descr)

	return (res)
}

id <- VAL_GSE_ID
betas <- loadSavedBetas(id)

rf_model <- loadSavedModel(REF_GSE_ID)
calfit <- loadSavedCalfit(REF_GSE_ID)

scores <- do.predict(rf_model, betas)

scores <- predict(
		calfit$glmnet.fit,
		newx = scores,
		type = "response",
		s = calfit$lambda.1se
	)[,, 1] # use lambda estimated by 10fold CVlambda

#res <- do.evaluateDIAG(scores)
res <- do.evaluateGEO(id, scores)

sum(res$MC == res$MC_PRED) / nrow(res)
sum(res$MCF == res$MCF_PRED) / nrow(res)


#write.xlsx(res, file="diag2_results.xlsx", rowNames=TRUE)

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



