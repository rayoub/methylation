
suppressWarnings(suppressPackageStartupMessages({
	library(conumee2)
	library(sesame)
	library(here)
}))

cnvPlotBatch <- function(batch_id) {
	output_dir = here::here("output", batch_id, "cnv")
	sample_dir <- here::here("data", "lab", batch_id)

	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    cnv_ref <- readRDS(file=here::here("results", "other", "CNV_REF_EPICv2.rds"))
    cnv_anno <- readRDS(file=here::here("results", "other", "CNV_ANNO_EPICv2.rds"))

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
			file.path(output_dir, paste0("cnv_", sample_id, ".pdf")),
			width = 12,
			height = 6
		)
		CNV.genomeplot(x, main = sample_id)
		dev.off()
	}
}