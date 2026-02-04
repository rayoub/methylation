
library(ranger)
library(openxlsx)
library(glmnet)
library(here)

source(here("R","constants.R"))
source(here("R","loading.R"))
source(here("R","mcf.R"))
source(here("R","mc.R"))

predictSampleScores <- function(betas) {

	# load rf model from reference data
	rf_model <- loadSavedModel(REF_GSE_ID)

	# predict raw scores from rf model
	rf_vars <- rf_model$forest$independent.variable.names
	filtered_betas <- betas[match(rf_vars, colnames(betas))]
	p <- predict(rf_model, data = filtered_betas)
	scores <- p$predictions
	rownames(scores) <- sub("_.*", "", rownames(filtered_betas))

	# calibrate scores using calfit model
	calfit <- loadSavedCalfit(REF_GSE_ID)
	scores <- predict(
			calfit$glmnet.fit,
			newx = scores,
			type = "response",
			s = calfit$lambda.1se
		)[,, 1] # use lambda estimated by 10fold CVlambda

	return(scores)
}

evaluateGEOSampleScores <- function (gse_id, scores) {

	anno <- loadSavedAnno(gse_id)
	if ("methylation class:ch1" %in% names(anno)) {

		res <- data.frame(MC = anno[,"methylation class:ch1"])

		mc_pred <- colnames(scores)[max.col(scores)]
		mcf <- unlist(ifelse(res[, "MC"] %in% names(MCF_LOOKUP), MCF_LOOKUP[res[, "MC"]], res[, "MC"]))
		res <- cbind(res, MC_PRED = mc_pred, MCF = mcf)

		mcf_pred <- unlist(ifelse(res[, "MC_PRED"] %in% names(MCF_LOOKUP), MCF_LOOKUP[res[, "MC_PRED"]], res[, "MC_PRED"]))
		res <- cbind(res, MCF_PRED = mcf_pred)

		mc_descr <- ifelse(res[, "MC"] %in% names(MC), MC[res[, "MC"]], NA)
		res <- cbind(res, MC_DESCR = mc_descr)
	}
	else {

		samples <- rownames(anno)	
		mc_pred <- colnames(scores)[max.col(scores)]
		res <- data.frame(SAMPLE_ID = samples, MC_PRED = mc_pred)

		mcf_pred <- unlist(ifelse(res[, "MC_PRED"] %in% names(MCF_LOOKUP), MCF_LOOKUP[res[, "MC_PRED"]], res[, "MC_PRED"]))
		res <- cbind(res, MCF_PRED = mcf_pred)

		mc_descr <- ifelse(res[, "MC_PRED"] %in% names(MC), MC[res[, "MC_PRED"]], NA)
		res <- cbind(res, MC_DESCR = mc_descr)

		nm_descr <- anno$`tumor type/grade:ch1`
		res <- cbind(res, NM_DESCR = nm_descr)
	}

	return (res)
}

evaluateDiagnosticSampleScores <- function (scores) {

	samples <- rownames(scores)	
	mc_pred <- colnames(scores)[max.col(scores)]
	res <- data.frame(SAMPLE_ID = samples, MC_PRED = mc_pred)

	mcf_pred <- unlist(ifelse(res[, "MC_PRED"] %in% names(MCF_LOOKUP), MCF_LOOKUP[res[, "MC_PRED"]], res[, "MC_PRED"]))
	res <- cbind(res, MCF_PRED = mcf_pred)

	mc_descr <- ifelse(res[, "MC_PRED"] %in% names(MC), MC[res[, "MC_PRED"]], NA)
	res <- cbind(res, MC_DESCR = mc_descr)

	return (res)
}
