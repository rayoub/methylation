# =============================
# 03 Train w/ Random Forest Classification
# =============================

library(randomForest)
library(glmnet)
library(matrixStats)

# beta_top: samples x probes (already top 32k variance probes) check
dim(beta_top)

# --- Step 1: Alignment ---
beta_gsm <- sub("_.*$", "", rownames(beta_top))   # Strip suffix if present
rownames(beta_top) <- beta_gsm

labels_filtered <- labels_filtered[labels_filtered$Sample %in% beta_gsm, ]
labels_filtered <- labels_filtered[match(rownames(beta_top), labels_filtered$Sample), ]

stopifnot(all(rownames(beta_top) == labels_filtered$Sample))
y <- as.factor(labels_filtered$Class)

# --- Step 2: Random Forest Feature Importance ---
n_probes <- ncol(beta_top)
chunk_size <- 10000
chunks <- split(1:n_probes, ceiling(seq_along(1:n_probes) / chunk_size))

importance_list <- list()

for (i in seq_along(chunks)) {
  message("Running chunk ", i, " of ", length(chunks), " (", length(chunks[[i]]), " probes)")
  
  rf <- randomForest(
    x = beta_top[, chunks[[i]], drop = FALSE],
    y = y,
    ntree = 100,
    mtry = 100,                        # same as Capper & NM (fixed at 100)
    sampsize = rep(8, length(levels(y))),  # downsample to 8 per class
    importance = TRUE
  )
  
  imp <- importance(rf, type = 1)      # Mean decrease in accuracy
  importance_list[[i]] <- imp
}

importance_all <- do.call(rbind, importance_list)

# --- Step 3: Rank and Select Top 10,000 Probes ---
# Rank each probe by its importance across all classes
imp_ranks <- apply(-importance_all, 2, rank, ties.method = "min")
min_rank <- apply(imp_ranks, 1, min)

top10k_idx <- order(min_rank)[1:10000]
beta_sel <- beta_top[, top10k_idx]

# --- Step 4: Final Random Forest Classifier ---
message("Training final random forest with top 10,000 probes...")

rf_final <- randomForest(
  x = as.data.frame(beta_sel),   # already samples × probes
  y = y,
  ntree = 10000,
  mtry = 100,
  sampsize = rep(8, length(levels(y))),
  importance = TRUE
)

# Force attach xnames if randomForest didn't store them
rf_final$xnames <- colnames(beta_sel)

saveRDS(rf_final, "rf_final_NMstyle2.rds")
saveRDS(top10k_idx, "rf_features_top10k.rds")



