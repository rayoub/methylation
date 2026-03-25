
suppressPackageStartupMessages({
	library(openxlsx)
	library(here)
})

source(here::here("R","files.R"))
source(here::here("R","preprocessing.R"))
source(here::here("R","mgmt.R"))
source(here::here("R","cnv.R"))
source(here::here("R","umap.R"))
source(here::here("R","prediction.R"))

processLabSamples <- function (batch_id) {

	#output_dir = here::here("output", batch_id)
	#dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

	#message("preprocessing samples ... ", Sys.time())
	#preprocessLabSamples(batch_id, "FFPE")

	message("performing classification ... ", Sys.time())
	batch_id <- "BATCH1"
	scores <- predictLabSampleScores(batch_id)
	#saveLabData(batch_id, "scores", scores)

	res <- evaluateLabSampleScores(scores)
	write.xlsx(res, file=here::here("output", batch_id, paste0(batch_id, ".xlsx")), rowNames=TRUE)

	#message("creating MGMT plots ... ", Sys.time())	
	#mgmtPlotBatch(batch_id)

	#message("creating CNV plots ... ", Sys.time())
	#cnvPlotBatch(batch_id)

	#message("creating UMAP plots ... ", Sys.time())
	#umapPlotBatch(batch_id)
	
	#message("processing samples complete ... ", Sys.time())
}