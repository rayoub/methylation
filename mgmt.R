# ==========================================
# MGMT Promoter Methylation Predictor
# ==========================================

# MGMTstp27 is not available on CRAN, download and install manually from source https://github.com/badozor/mgmtstp27
library(mgmtstp27)
library(minfi)
library(dplyr)
library(ggtext)

mgmt.plot <- function(diag_id, sample_id, preds) {
  cutoff <- 0.3582 # STP27 cutoff - This is a hardcoded value biologically
  x_max <- 1.5 # extra space for right-hand annotation

  # --- Extract sample data ---
  df_one <- preds %>%
    filter(sample == sample_id) %>%
    mutate(
      Class = ifelse(state == "M", "Methylated", "Unmethylated")
    )

  # --- Styled annotation text ---
  info_text_html <- paste0(
    "<b>Status:</b> ",
    df_one$Class,
    "<br>",
    "<b>Estimate:</b> ",
    sprintf("%.6f", df_one$pred),
    "<br>",
    "<b>CI-Lower:</b> ",
    sprintf("%.6f", df_one$lower),
    "<br>",
    "<b>CI-Upper:</b> ",
    sprintf("%.6f", df_one$upper),
    "<br>",
    "<b>CI-Cutoff:</b> ",
    sprintf("%.6f", cutoff)
  )

  # --- Plot ---
  ggplot(df_one, aes(y = sample)) +
    # Background first (so white lines appear on top)
    geom_rect(
      aes(xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf),
      fill = "gray95",
      color = NA
    ) +
    # White vertical guide lines (now visible)
    geom_vline(
      xintercept = c(0, 0.25, 0.5, 0.75, 1),
      color = "white",
      linewidth = 0.72,
      alpha = 0.8
    ) +
    # Confidence interval
    geom_errorbar(
      aes(xmin = lower, xmax = upper),
      width = 0.1,
      color = "gray40",
      linewidth = 1,
      orientation = "y"
    ) +
    # Prediction point
    geom_point(aes(x = pred, color = Class), size = 5) +
    # Cutoff line
    geom_vline(
      xintercept = cutoff,
      linetype = "solid",
      color = "red",
      linewidth = 1
    ) +
    # Axis range + clipping
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    # Axis ticks
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1.0")
    ) +
    # Colors
    scale_color_manual(
      values = c("Methylated" = "#D73027", "Unmethylated" = "#4575B4")
    ) +
    # Theme
    theme_minimal(base_size = 12) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = ggplot2::margin(
        t = 10,
        r = 140,
        b = 10,
        l = 10,
        unit = "pt"
      ),
      axis.ticks.x = element_line(color = "white", linewidth = 0.5),
      axis.ticks.length.x = unit(3, "pt")
    ) +
    labs(
      title = "MGMT Promoter Status Prediction (MGMT-STP27)",
      x = sprintf("Score (red line cutoff = %.4f)", cutoff)
    ) +
    # Right-side annotation
    annotate(
      "richtext",
      x = 1.02,
      y = 1,
      label = info_text_html,
      hjust = 0,
      vjust = 0.5,
      size = 4.4,
      fill = NA,
      label.color = NA,
      family = "sans",
      lineheight = 1.2
    )

    ggsave(file.path("output", diag_id, paste0("MGMT_", sample_id, ".pdf")))
}

# gather input
diag_id <- "DIAG1"

# output directory
dir.create(file.path("output", diag_id), showWarnings = FALSE)

# read meth arrays
data_dir <- file.path("data", "diagnostic", diag_id)
rgset <- read.metharray.exp(data_dir, verbose = TRUE)
mset <- preprocessRaw(rgset)

# convert M values
mvals <- log2((getMeth(mset) + 1) / (getUnmeth(mset) + 1))
mvals <- as.data.frame(t(mvals))

# predict 
preds <- MGMTpredict(mvals)

# plot
for (sample_id in rownames(preds)) {
  mgmt.plot(diag_id, sample_id, preds)
}


