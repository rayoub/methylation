
suppressWarnings(suppressPackageStartupMessages({
	library(ranger)
	library(glmnet)
	library(tidyr)
	library(here)
}))

source(here::here("R","constants.R"))
source(here::here("R","files.R"))
source(here::here("R","mcf.R"))
source(here::here("R","mc.R"))

predictSampleScores <- function(betas) {

	# load rf model from reference data
	rf_model <- loadGeoData(REF_GSE_ID, "model")

	# predict raw scores from rf model
	rf_vars <- rf_model$forest$independent.variable.names
	filtered_betas <- betas[match(rf_vars, colnames(betas))]
	p <- predict(rf_model, data = filtered_betas)
	scores <- p$predictions
	rownames(scores) <- rownames(filtered_betas)

	# calibrate scores using calfit model
	calfit <- loadGeoData(REF_GSE_ID, "calfit")
	scores <- predict(
			calfit$glmnet.fit,
			newx = scores,
			type = "response",
			s = calfit$lambda.1se
		)[,, 1] # use lambda estimated by 10fold CVlambda

	return(scores)
}

predictGeoSampleScores <- function (gse_id) {
	
	betas <- loadGeoData(gse_id, "betas")
	predictSampleScores(betas)
}

predictLabSampleScores <- function (batch_id) {
	
	betas <- loadLabData(batch_id, "betas")
	predictSampleScores(betas)
}

evaluateGEOSampleScores <- function (scores, gse_id) {

	anno <- loadGeoData(gse_id, "anno")
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

evaluateLabSampleScores <- function (scores) {

	scores_df <- as.data.frame(scores, row.names = rownames(scores))
	scores_t <- tibble::rownames_to_column(scores_df, var = "sample_id")
	
	scores_mc <- scores_t |>
		pivot_longer(
			cols = !(sample_id),
			names_to = "mc",
			values_to = "mc_score"
		) |>
		dplyr::group_by(sample_id) |>
		dplyr::slice_max(mc_score) 

	scores_mcf <- scores_t |>
		pivot_longer(
			cols = !(sample_id),
			names_to = "mc",
			values_to = "mc_score"
		) |>
		dplyr::mutate(
			mcf = mcf_lookup(mc),
			.after = mc
		) |>
		dplyr::group_by(sample_id, mcf) |>
		dplyr::summarize(
			mcf_score = sum(mc_score),
			.groups = "drop_last"
		) |>
		dplyr::slice_max(order_by = mcf_score)

	scores_j <- scores_mc |> 
		dplyr::inner_join(scores_mcf, join_by(sample_id)) |>
		dplyr::mutate(
			mc_descr = MC[mc]
		)

	return (scores_j)
}
