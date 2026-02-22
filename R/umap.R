suppressWarnings(suppressPackageStartupMessages({
  library(uwot)
  library(ggplot2)
  library(GEOquery)
  library(dplyr)
  library(ggrepel)
  library(here)
}))

source(here::here("R", "constants.R"))
source(here::here("R", "files.R"))

umapPlot <- function(umap_betas, umap_model, umap_df, centroids, sample_id, sample) {
  # calculate min and max values
  x_min <- min(umap_df$UMAP1)
  x_max <- max(umap_df$UMAP1)
  y_min <- min(umap_df$UMAP2)
  y_max <- max(umap_df$UMAP2)

  x_cut = x_min + 0.75 * (x_max - x_min)
  y_cut = y_min + 0.75 * (y_max - y_min)

  # umap sample
  coords <- umap_transform(sample, model = umap_model)
  umap_sample <- as.data.frame(coords)
  colnames(umap_sample) <- c("UMAP1", "UMAP2")
  umap_sample$Class <- sample_id

  # directions
  x <- umap_sample$UMAP1
  y <- umap_sample$UMAP2
  x_sign <- if_else(x < x_cut, 1, -1)
  y_sign <- if_else(y < y_cut, 1, -1)

  # plot
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Class)) +
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
      data = umap_sample,
      aes(x = UMAP1 + (x_sign * 3), y = UMAP2 + (y_sign * 3), label = Class),
      size = 4,
      fontface = "bold",
      color = "black",
      show.legend = FALSE
    ) +
    geom_segment(
      data = umap_sample,
      aes(
        x = UMAP1,
        y = UMAP2,
        xend = UMAP1 + (x_sign * 2.4),
        yend = UMAP2 + (y_sign * 2.4)
      ),
      linewidth = 0.5,
      color = "black",
      arrow = arrow(
        length = unit(0.1, "cm"),
        type = "closed",
        ends = "first"
      ),
      show.legend = FALSE
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = ggplot2::margin(
        t = 10,
        r = 10,
        b = 10,
        l = 10,
        unit = "pt"
      )
    ) +
    labs(title = "GSE90496 Reference Set")
}

umapPlotBatch <- function(batch_id) {
  # output directory
  output_dir = here::here("output", batch_id, "umap")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # load umap_betas, umap_df, centroids
  load(here::here("results", "geo", paste0(REF_GSE_ID, "_umap.RData")))

  # load umap for ref data
  umap_model <- umap(
    umap_betas,
    n_neighbors = 30,
    min_dist = 0.75,
    metric = "euclidean",
    ret_model = TRUE
  )

  # load batch data
  betas <- loadLabData(batch_id, "betas")

  # iterate samples
  for (sample_id in rownames(betas)) {
    # set seed
    seed <- 180314
    set.seed(seed, kind = "L'Ecuyer-CMRG")

    sample <- betas[sample_id, colnames(umap_betas)]

    p <- umapPlot(umap_betas, umap_model, umap_df, centroids, sample_id, sample)
    ggsave(
      file.path(output_dir, paste0("umap_", sample_id, ".pdf"))
    )
  }
}
