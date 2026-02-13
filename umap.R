
library(uwot)
library(ggplot2)
library(GEOquery)
library(dplyr)
library(ggrepel)
library(here)

source(here("R", "constants.R"))
source(here("R", "loading.R"))

# set seed
seed <- 180314
set.seed(seed, kind = "L'Ecuyer-CMRG")

# load ref data
umap_betas <- loadSavedUmapBetas(REF_GSE_ID)
anno <- loadSavedAnno(REF_GSE_ID)
labels <- data.frame(
  Sample = anno$geo_accession,
  Class = anno$`methylation class:ch1`
)

# load diag data
diag_id <- "DIAG1"
betas <- loadSavedBetas(diag_id)

# iterate samples
sample_id <- "209794430029_R05C01"
sample <- betas[sample_id,colnames(umap_betas)]

# umap
augmented_betas <- rbind(umap_betas, sample)
umap_res <- umap(
  augmented_betas,
  n_neighbors = 30,
  min_dist = 0.75,
  metric = "euclidean"
)

# set names 
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_first <- umap_df[-nrow(umap_df),]
rownames(umap_first) <- labels$Sample
umap_first$Class <- labels$Class

# calculate centroids for labels
centroids <- umap_first |>
  group_by(Class) |>
  summarize(
    UMAP1 = mean(UMAP1),
    UMAP2 = mean(UMAP2)
  )

umap_last <- umap_df[nrow(umap_df),]
umap_last$Class = sample_id

# plot
ggplot(umap_first, aes(x = UMAP1, y = UMAP2, color = Class)) +
  geom_point(
    size = 0.9, 
    alpha = 0.5, 
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = centroids,
    aes(label = Class),
    size = 2,
    fontface = "bold",
    color = "black",
    max.overlaps = Inf
  ) +
  geom_label(
    data = umap_last, 
    aes(x = UMAP1 + 3, y = UMAP2 + 3, label = Class),
    size = 4,
    fontface = "bold",
    color = "black",
    show.legend = FALSE
  ) +
  geom_segment(
    data = umap_last, 
    aes(x = UMAP1, y = UMAP2, xend = UMAP1 + 2.5, yend = UMAP2 + 2.5),
    linewidth = 1,
    color = "black",
    arrow = arrow(length = unit(0.2, "cm"), type = "closed", ends = "first"),
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(title = "GSE90496 Reference Set")