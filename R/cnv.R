
suppressWarnings(suppressPackageStartupMessages({
	library(conumee2)
	library(sesame)
	library(here)
}))

source(here::here("R", "loading.R"))

cnvPlotBatch <- function(diag_id) {
	output_dir = here::here("output", diag_id, "cnv")
	sample_dir <- here::here("data", diag_id)

	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

	cnv_ref <- loadSavedCnvRef()
	cnv_anno <- loadSavedCnvAnno()

	message("open sesame ... ", Sys.time())
	cnv_query <- openSesame(sample_dir, prep = "QCDPB", func = NULL)

	message("genome plots ... ", Sys.time())
	for (sample_id in names(cnv_query)) {
		cnv_data <- CNV.load(totalIntensities(cnv_query[[sample_id]]))

		x <- CNV.fit(cnv_data, cnv_ref, cnv_anno)
		x <- CNV.bin(x)
		x <- CNV.detail(x)
		x <- CNV.segment(x)

		pdf(
			file.path(output_dir, paste0(sample_id, ".pdf")),
			width = 12,
			height = 6
		)
		CNV.genomeplot(x, main = sample_id)
		dev.off()
	}
}
