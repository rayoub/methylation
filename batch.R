
suppressPackageStartupMessages({
	library(openxlsx)
	library(here)
})

source(here::here("R","loading.R"))
source(here::here("R","preprocessing.R"))
source(here::here("R","mgmt.R"))
source(here::here("R","cnv.R"))
source(here::here("R","umap.R"))
source(here::here("R","prediction.R"))

processSamples <- function (diag_id) {

	output_dir = here::here("output", diag_id)
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

	message("preprocessing samples ... ", Sys.time())
	preprocessDiagnosticSamples(diag_id, "FFPE")

	message("performing classification ... ", Sys.time())
	scores <- predictSampleScores(diag_id)
	res <- evaluateDiagnosticSampleScores(scores)
	write.xlsx(res, file=here::here("output", diag_id, paste0(diag_id, ".xlsx")), rowNames=TRUE)

	message("creating MGMT plots ... ", Sys.time())	
	mgmtPlotBatch(diag_id)

	message("creating CNV plots ... ", Sys.time())
	cnvPlotBatch(diag_id)

	message("creating UMAP plots ... ", Sys.time())
	umapPlotBatch(diag_id)
	
	message("processing samples complete ... ", Sys.time())
}