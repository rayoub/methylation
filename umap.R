
library(uwot)
library(ggplot2)
library(GEOquery)
library(dplyr)
library(ggrepel)
library(here)

source(here("R", "constants.R"))
source(here("R", "loading.R"))

N_TOP_BETAS <- 32000 # top variably methylated probes

# load betas and filter top 32,000 variably methylated probes
betas <- loadSavedBetas(REF_GSE_ID)
betas <- betas[, order(-apply(betas, 2, sd))[1:N_TOP_BETAS]]

# load sample names and classes
anno <- loadSavedAnno(REF_GSE_ID)
labels <- data.frame(
  Sample = anno$geo_accession,
  Class = anno$`methylation class:ch1`
)

# set seed
seed <- 180314
set.seed(seed, kind = "L'Ecuyer-CMRG")

# umap
umap_res <- umap(
  betas,
  n_neighbors = 30,
  min_dist = 0.75,
  metric = "euclidean"
)

umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
rownames(umap_df) <- labels$Sample
umap_df$Class <- labels$Class

centroids <- umap_df |>
  group_by(Class) |>
  summarize(
    UMAP1 = mean(UMAP1),
    UMAP2 = mean(UMAP2)
  )

# plot
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Class)) +
  geom_point(size = 0.9, alpha = 0.5) +
  geom_text_repel(
    data = centroids,
    aes(label = Class),
    size = 2,
    fontface = "bold",
    color = "black",
    max.overlaps = Inf
  ) +
  guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))) +
  theme_minimal() +
  labs(title = "GSE90496 Reference Set") +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.box = "horizontal"
  )
