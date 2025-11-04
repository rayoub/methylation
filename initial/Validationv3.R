# ---- Libraries ----
library(GEOquery)
library(minfi)
library(preprocessCore)
library(randomForest)
library(glmnet)
library(matrixStats)
library(R.utils)
library(dplyr)
library(stringr)
library(readr)

set.seed(42)

# ---- Parameters ----
gse_val_id <- "GSE109379"
n_val_samples <- 50       # how many GSMs to process

# ---- 0) Load training artifacts ----
beta_train <- readRDS("data_raw/beta_filtered_GSE90496.RDS")
top10k_idx <- readRDS("rf_features_top10k.rds")
rf_final <- readRDS("rf_final_full2.rds")
glmnet_final <- readRDS("glmnet_recalibration_NMstyle.rds")

# Rebuild top32k and top10k probe IDs
vars_all <- matrixStats::rowVars(beta_train, na.rm = TRUE)
top32k_idx <- order(vars_all, decreasing = TRUE)[1:32000]
beta_top32k <- beta_train[top32k_idx, , drop = FALSE]
top10k_probes <- rownames(beta_top32k)[top10k_idx]
cat("Top10k probes reconstructed:", length(top10k_probes), "\n")

# ---- 1) Download + preprocess validation data ----
gse_val <- getGEO(gse_val_id, GSEMatrix = FALSE)
gsm_list_val <- GSMList(gse_val)
val_samples <- names(gsm_list_val)[1:n_val_samples]

# Collect supplementary URLs
supp_urls_val <- list()
for (gsm_name in val_samples) {
  sfiles <- Meta(gsm_list_val[[gsm_name]])$supplementary_file
  if (length(sfiles) > 0) supp_urls_val[[gsm_name]] <- sfiles
}
unique_supp_val <- gsub("^ftp://", "https://", unique(unlist(supp_urls_val)))

dir.create("data_raw/validation_idats", showWarnings = FALSE, recursive = TRUE)
for (url in unique_supp_val) {
  fname <- file.path("data_raw/validation_idats", basename(url))
  if (!file.exists(fname)) {
    message("Downloading ", basename(url))
    try(download.file(url, destfile = fname, mode = "wb", method = "curl"))
    if (grepl("\\.zip$", fname)) unzip(fname, exdir = "data_raw/validation_idats")
    if (grepl("\\.tar|\\.tgz$", fname)) untar(fname, exdir = "data_raw/validation_idats")
  }
}

# Unzip .gz files
gz_files <- list.files("data_raw/validation_idats", pattern = "idat.gz$", full.names = TRUE, recursive = TRUE)
for (f in gz_files) R.utils::gunzip(f, overwrite = TRUE)

# ---- 2) Preprocess with Funnorm ----
rgSet <- read.metharray.exp("data_raw/validation_idats", recursive = TRUE)
mset <- preprocessFunnorm(rgSet)
beta_val_full <- getBeta(mset)
cat("Validation beta dims:", dim(beta_val_full), "\n")

# ---- 3) Subset to top10k probes, impute, and QN ----
present <- top10k_probes %in% rownames(beta_val_full)
cat("Top10k present in validation:", sum(present), "of", length(top10k_probes), "\n")

beta_val_top10k <- matrix(NA, nrow = length(top10k_probes), ncol = ncol(beta_val_full),
                          dimnames = list(top10k_probes, colnames(beta_val_full)))
have <- intersect(top10k_probes, rownames(beta_val_full))
beta_val_top10k[have, ] <- beta_val_full[have, , drop = FALSE]

# Impute missing probes with training means
beta_train_top10k <- beta_top32k[top10k_idx, , drop = FALSE]
tr_means <- rowMeans(beta_train_top10k, na.rm = TRUE)
missing_probes <- setdiff(top10k_probes, have)
if (length(missing_probes) > 0) {
  message("Imputing ", length(missing_probes), " missing probes with training probe means.")
  beta_val_top10k[missing_probes, ] <- matrix(tr_means[missing_probes],
                                              nrow = length(missing_probes),
                                              ncol = ncol(beta_val_full), byrow = FALSE)
}

# Quantile normalize (validation -> training distribution)
cat("Constructing quantile target from training top10k...\n")
sorted_train_cols <- apply(beta_train_top10k, 2, sort, na.last = NA)
train_ref <- rowMeans(sorted_train_cols, na.rm = TRUE)
beta_val_qn <- normalize.quantiles.use.target(as.matrix(beta_val_top10k), target = train_ref)
rownames(beta_val_qn) <- rownames(beta_val_top10k)
colnames(beta_val_qn) <- colnames(beta_val_top10k)
cat("QN complete. val qn range:", range(beta_val_qn, na.rm = TRUE), "\n")

# ---- 5) Prepare RF input (beta values, no z-score) ----
beta_val_for_rf <- t(beta_val_final)
rownames(beta_val_for_rf) <- colnames(beta_val_final)
colnames(beta_val_for_rf) <- rownames(beta_val_final)

stopifnot(all(colnames(beta_val_for_rf) == rf_final$xnames))
cat("RF input dims:", dim(beta_val_for_rf), " Range:", range(beta_val_for_rf), "\n")

# ---- 6) Predict RF ----
rf_probs_val <- predict(rf_final, newdata = as.data.frame(beta_val_for_rf), type = "prob")
cat("RF prob matrix dims:", dim(rf_probs_val), "\n")
print(summary(as.vector(rf_probs_val)))


# ---- 7) Calibration (glmnet) ----
if (!is.null(glmnet_final)) {
  lambda_to_use <- if (!is.null(glmnet_final$lambda.1se)) glmnet_final$lambda.1se else glmnet_final$lambda.min
  rf_probs_mat <- as.matrix(rf_probs_val)
  calib_array <- predict(glmnet_final, newx = rf_probs_mat, s = lambda_to_use, type = "response")
  if (is.array(calib_array)) {
    calibrated_probs_val <- calib_array[,,1]
  } else if (is.list(calib_array) && length(calib_array) > 0) {
    calibrated_probs_val <- do.call(cbind, lapply(calib_array, function(x) x[,1]))
  } else {
    calibrated_probs_val <- calib_array
  }
  colnames(calibrated_probs_val) <- colnames(rf_probs_val)
} else {
  calibrated_probs_val <- as.matrix(rf_probs_val)
}

# Using raw RF instead of calibration to troubleshoot - I think the calibration step is messed up
calibrated_probs_val <- rf_probs_val



# ---- 8) Results ----
pred_class_val <- colnames(calibrated_probs_val)[apply(calibrated_probs_val, 1, which.max)]
pred_class_prob_val <- apply(calibrated_probs_val, 1, max)

results_val <- data.frame(
  Sample = rownames(calibrated_probs_val),
  Predicted_Class = pred_class_val,
  Calibrated_Probability = pred_class_prob_val,
  stringsAsFactors = FALSE
)




cat("\nClassification confidence summary:\n")
print(table(results_val$Confidence))
cat("Results (first few):\n"); print(head(results_val))
print(results_val)
cat("Summary of calibrated probabilities:\n"); print(summary(results_val$Calibrated_Probability))



