# =============================
# 05 Evaluate Model Calibration
# =============================

library(randomForest)
library(glmnet)
library(pROC)

# Predicted class = highest probability
pred_class <- colnames(rf_preds)[max.col(rf_preds, ties.method = "first")]

conf_mat <- caret::confusionMatrix(
  factor(pred_class, levels = levels(y)),
  y
)
print(conf_mat)

# ROC / log-loss
roc_list <- multiclass.roc(y, rf_preds)
print(roc_list$auc)
#This replicates the same AUC from the Capper paper


# Now train on all 2800 samples:  Final RF on all data
rf_final_full <- randomForest(
  x = beta_sel,
  y = y,
  ntree = 10000,
  mtry = 100,
  sampsize = rep(8, length(levels(y)))
)
saveRDS(rf_final_full, "rf_final_full2.rds")

# Predict RF probabilities for all samples
rf_scores_all <- predict(rf_final_full, beta_sel, type = "prob")

# Fit logistic regression recalibration on full data
glmnet_final <- glmnet(
  x = rf_scores_all,
  y = y,
  family = "multinomial",
  alpha = 0
)
saveRDS(glmnet_final, "glmnet_final_full2.rds")
