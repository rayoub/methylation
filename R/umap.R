suppressWarnings(suppressPackageStartupMessages({
  library(uwot)
  library(ggplot2)
  library(GEOquery)
  library(dplyr)
  library(ggrepel)
  library(here)
}))

source(here("R", "constants.R"))
source(here("R", "loading.R"))

umapPlotBatch <- function(diag_id) {

  # output directory
  output_dir = here("output", diag_id, "umap")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # load ref data
  umap_betas <- loadSavedUmapBetas(REF_GSE_ID)
  anno <- loadSavedAnno(REF_GSE_ID)
  labels <- data.frame(
    Sample = anno$geo_accession,
    Class = anno$`methylation class:ch1`
  )

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

  # calculate min and max values
  x_min <- min(umap_df$UMAP1)
  x_max <- max(umap_df$UMAP1)
  y_min <- min(umap_df$UMAP2)
  y_max <- max(umap_df$UMAP2)

  x_cut = x_min + 0.75 * (x_max - x_min)
  y_cut = y_min + 0.75 * (y_max - y_min)
    
  # calculate ref data umap centroids 
  centroids <- umap_df |>
    group_by(Class) |>
    summarize(
      UMAP1 = mean(UMAP1),
      UMAP2 = mean(UMAP2)
    )

  # load diag data
  betas <- loadSavedBetas(diag_id)

  # iterate samples
  for (sample_id in rownames(betas)) {
  
    # set seed
    seed <- 180314
    set.seed(seed, kind = "L'Ecuyer-CMRG")

    sample <- betas[sample_id, colnames(umap_betas)]

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
        aes(x = UMAP1, y = UMAP2, xend = UMAP1 + (x_sign * 2.4), yend = UMAP2 + (y_sign * 2.4)),
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

    ggsave(
      file.path(output_dir, paste0("UMAP_", sample_id, ".pdf"))
    )
  }
}
