# =============================
# 04 Score Recalibration
# =============================
# Turn random forest raw scores into calibrated class probabilities

library(randomForest)
library(glmnet)
library(caret)

set.seed(42)

# Load final feature subset (10k) and labels
top10k_idx <- readRDS("rf_features_top10k.rds")
beta_sel <- beta_top[, top10k_idx]
y <- as.factor(labels_filtered$Class)

# 5-fold CV
n_folds <- 5
folds <- createFolds(y, k = n_folds, list = TRUE, returnTrain = TRUE)

# Placeholder for out-of-fold predictions
rf_preds <- matrix(NA, nrow = length(y), ncol = length(levels(y)))
colnames(rf_preds) <- levels(y)

for (i in seq_along(folds)) {
  message("Fold ", i, "/", n_folds)
  
  train_idx <- folds[[i]]
  test_idx <- setdiff(seq_along(y), train_idx)
  
  # --- Adjust sampsize per class to avoid error ---
  class_counts <- table(y[train_idx])
  sampsize_adj <- pmin(class_counts, 8)  # 8 or fewer if class has <8 samples
  
  rf <- randomForest(
    x = beta_sel[train_idx, , drop = FALSE],
    y = y[train_idx],
    ntree = 10000,
    mtry = 100,
    sampsize = sampsize_adj
  )
  
  # Predict probabilities for held-out fold
  fold_pred <- predict(rf, beta_sel[test_idx, , drop = FALSE], type = "prob")
  rf_preds[test_idx, ] <- fold_pred
}

# --- Train L2-penalized multinomial logistic regression ---
x <- rf_preds
y_glm <- y

glmnet_model <- cv.glmnet(
  x = x,
  y = y_glm,
  family = "multinomial",
  alpha = 0,       # L2 penalty
  nfolds = 5,
  type.measure = "deviance"
)

saveRDS(glmnet_model, "glmnet_recalibration_NMstyle2.rds")
message("✅ Recalibration model trained and saved.")