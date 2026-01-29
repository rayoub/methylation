## -- 02 Visualize Dimensional reduction -- ##

# Load librarires
library(matrixStats)
library(uwot)
library(ggplot2)
library(cluster)
library(ggforce)
library(RColorBrewer)
library(GEOquery)
library(viridis)
library(Polychrome)
library(dplyr)
library(ggrepel)

# -----------------------------
# 1. Load beta matrix
# -----------------------------
beta <- readRDS("data_raw/beta_filtered_GSE90496.RDS")
dim(beta)

# -----------------------------
# 2. Select top variable probes
# -----------------------------
vars <- rowVars(beta, na.rm = TRUE)
top_idx <- order(vars, decreasing = TRUE)[1:32000]   # top 32k probes like Capper
beta_top <- t(beta[top_idx, ])   # transpose → samples x probes

# -----------------------------
# 3. Load metadata
# -----------------------------
gse_list <- getGEO("GSE90496", GSEMatrix = TRUE)
gse <- if(is.list(gse_list)) gse_list[[1]] else gse_list
meta <- pData(gse)

# Extract sample names using geo_accession
labels <- data.frame(
  Sample = meta$geo_accession,
  Class = meta$`methylation class:ch1`
)

# -----------------------------
# 4. Match beta matrix columns to metadata GSM IDs
# -----------------------------
# Extract GSM IDs from beta column names (strip chip suffix)
beta_gsm <- sub("_.*$", "", rownames(beta_top))

# Find common samples
common_samples <- intersect(beta_gsm, labels$Sample)
length(common_samples)  # should be 2801

# Subset beta_top and labels to common samples
beta_top <- beta_top[beta_gsm %in% common_samples, ]
labels_filtered <- labels[labels$Sample %in% common_samples, ]
labels_filtered <- labels_filtered[match(beta_gsm[beta_gsm %in% common_samples], labels_filtered$Sample), ]

# -----------------------------
# 5. Run UMAP
# -----------------------------
set.seed(123)
n_samples <- nrow(beta_top)
n_neighbors <- min(30, n_samples - 1)  # safe default

umap_res <- umap(
  beta_top,
  n_neighbors = n_neighbors,
  min_dist = 0.75,
  metric = "euclidean"
)

umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample <- rownames(beta_top)

# Merge with labels
umap_df$Class <- labels_filtered$Class

num_classes <- length(unique(umap_df$Class))

centroids <- umap_df %>%
  group_by(Class) %>%
  summarize(
    UMAP1 = mean(UMAP1),
    UMAP2 = mean(UMAP2)
  )

# -----------------------------
# 6. Plot UMAP
# -----------------------------
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


