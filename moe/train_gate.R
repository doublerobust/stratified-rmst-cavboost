#!/usr/bin/env Rscript
# Phase 1: Train RF gate on full training set, save for chunk evaluation
#
# Usage:
#   Rscript moe/train_gate.R
#
# Output:
#   moe/results/rds/gate_model.rds      — trained ranger RF
#   moe/results/rds/gate_test_data.rds  — test set rows + metadata
#
# Then run:
#   python3 moe/run_parallel_gate.py    — runs Phase 2 in parallel chunks

library(ranger)
library(data.table)

source("moe/scenario_generator.R")

.cov <- function(d) setdiff(names(d), c("trt01p","time","status","A","U","delta_tilde"))

cat("Loading gate training data...\n")
d <- fread("moe/results/gate_training_data.csv")

exclude_cols <- c("seed","family","n_train","optimal_K","optimal_method","oracle_rate",
  "auc_K1","auc_K2","auc_K3","auc_K4","auc_K5","auc_VT",
  "cfg_n_predictive","cfg_n_prognostic","cfg_te_scale","cfg_b0",
  "cfg_prognostic_form","cfg_te_start","cfg_te_peak","cfg_te_decay",
  grep("^K[1-5]_", names(d), value = TRUE))
feat_cols <- setdiff(names(d), exclude_cols)
is_num <- sapply(d[, ..feat_cols], is.numeric)
feat_cols <- feat_cols[is_num]
cat(sprintf("  %d rows, %d features\n", nrow(d), length(feat_cols)))

X <- as.matrix(d[, ..feat_cols])
for (j in seq_len(ncol(X))) X[is.na(X[,j]), j] <- mean(X[,j], na.rm=TRUE)
X[!is.finite(X)] <- 0
y <- d$optimal_method

# ---- Train/test split at config level ----
set.seed(20260530)
if (length(unique(paste0(d$family, d$n_train))) < nrow(d)) {
  config_id <- rleidv(d[, c("family", "n_train")])
  unique_configs <- unique(config_id)
  n_configs <- length(unique_configs)
  train_configs <- sample(unique_configs, round(0.8 * n_configs))
  idx <- which(config_id %in% train_configs)
} else {
  idx <- sample(nrow(X), round(0.8 * nrow(X)))
}
X_tr <- X[idx,]; y_tr <- y[idx]
X_te <- X[-idx,]; y_te <- y[-idx]
d_test <- d[-idx, ]
cat(sprintf("  Train: %d, Test: %d\n", length(y_tr), length(y_te)))

# ---- Train Random Forest gate ----
cat("Training gate (ranger random forest)...\n")
gate <- ranger(
  x = X_tr, y = factor(y_tr),
  num.trees = 500,
  mtry = floor(sqrt(ncol(X_tr))),
  importance = "impurity",
  seed = 20260530,
  probability = TRUE
)

# Gate predictions on test set
probs <- predict(gate, X_te)$predictions
gate_K <- max.col(probs, ties.method = "first")

# ---- Save ----
dir.create("moe/results/rds", showWarnings = FALSE, recursive = TRUE)
saveRDS(gate, "moe/results/rds/gate_model.rds")
saveRDS(list(X_te = X_te, d_test = d_test, gate_K = gate_K),
        "moe/results/rds/gate_test_data.rds")
cat(sprintf("Saved: moe/results/rds/gate_model.rds + gate_test_data.rds\n"))
cat(sprintf("Train: %d rows, Test: %d reps across %d configs\n",
    nrow(X_tr), nrow(d_test), length(unique(paste0(d_test$family, d_test$n_train)))))
