
library(conumee2)
library(sesame)
library(here)

source(here("R", "loading.R"))

diag_id <- "DIAG1"

sample_dir <- file.path("data", "diagnostic", diag_id)

cnv_ref <- loadSavedCnvRef()
cnv_anno <- loadSavedCnvAnno() 

cnv_query <- openSesame(sample_dir, prep = "QCDPB", func = NULL)

for (sample_id in names(cnv_query)) {

	cnv_data <- CNV.load(totalIntensities(cnv_query[[sample_id]]))

	x <- CNV.fit(cnv_data, cnv_ref, cnv_anno)
	x <- CNV.bin(x)
	x <- CNV.detail(x)	
	x <- CNV.segment(x)

	pdf(file.path("output", paste0(sample_id, ".pdf")), width = 12, height = 6)
	CNV.genomeplot(x)
	dev.off()

	break
}
