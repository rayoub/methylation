library(uwot)
library(here)
library(dplyr)

source(here::here("R", "constants.R"))
source(here::here("R", "files.R"))

N_TOP_BETAS <- 32000 # top variably methylated probes

# load betas and annotation
betas <- loadGeoData(REF_GSE_ID, "betas")
anno <- loadGeoData(REF_GSE_ID, "anno")
labels <- data.frame(
	Sample = anno$geo_accession,
	Class = anno$`methylation class:ch1`
)

# filter top 32,000 variably methylated probes
umap_betas <- betas[, order(-apply(betas, 2, sd))[1:N_TOP_BETAS]]

# load umap for ref data
umap_model <- umap(
	umap_betas,
	n_neighbors = 30,
	min_dist = 0.75,
	metric = "euclidean",
	ret_model = TRUE
)

# set names
umap_df <- as.data.frame(umap_model$embedding)
colnames(umap_df) <- c("UMAP1", "UMAP2")
rownames(umap_df) <- labels$Sample
umap_df$Class <- labels$Class

# calculate ref data umap centroids
centroids <- umap_df |>
	group_by(Class) |>
	summarize(
		UMAP1 = mean(UMAP1),
		UMAP2 = mean(UMAP2)
	)

save(
	umap_betas,
	umap_df,
	centroids,
	file = here::here("results", paste0(REF_GSE_ID, "_umap.RData"))
)
