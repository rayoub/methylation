
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
filtered_betas <- loadSavedBetas(REF_GSE_ID)
filtered_betas <- filtered_betas[, order(-apply(filtered_betas, 2, sd))[1:N_TOP_BETAS]]

saveRDS(filtered_betas, file=here("results", paste0(REF_GSE_ID, "_filtered_betas.rds")))

anno <- loadSavedAnno(REF_GSE_ID)

diag_id <- "DIAG1"
betas <- loadSavedBetas(diag_id)

sample_id <- "209794430029_R05C01"
sample <- betas[sample_id,colnames(filtered_betas)]

# load sample names and classes
labels <- data.frame(
  Sample = c(anno$geo_accession, sample_id),
  Class = c(anno$`methylation class:ch1`, "CURRENT SAMPLE")
)

augmented_betas <- rbind(filtered_betas, sample)

# args to function (sample, labels)



# set seed
seed <- 180314
set.seed(seed, kind = "L'Ecuyer-CMRG")

# umap
umap_res <- umap(
  augmented_betas,
  n_neighbors = 30,
  min_dist = 0.75,
  metric = "euclidean"
)

# set names 
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
rownames(umap_df) <- labels$Sample
umap_df$Class <- labels$Class

# calculate centroids for labels
centroids <- umap_df |>
  group_by(Class) |>
  summarize(
    UMAP1 = mean(UMAP1),
    UMAP2 = mean(UMAP2)
  )



# plot
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Class)) +
  geom_point(size = 0.9, alpha = 0.5, show.legend = FALSE) +
  geom_text_repel(
    data = centroids,
    aes(label = Class),
    size = 2,
    fontface = "bold",
    color = "black",
    max.overlaps = Inf
  ) +
  #guides(color = guide_legend(ncol = 3, override.aes = list(size = 3))) +
  theme_minimal() +
  labs(title = "GSE90496 Reference Set") +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.box = "horizontal"
  )
