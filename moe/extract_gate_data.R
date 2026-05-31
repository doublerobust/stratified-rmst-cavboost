#!/usr/bin/env Rscript
# Extract gate training data from RDS results + raw data
# Dumps features + optimal K + richer meta-features to CSV
library(data.table)
library(glmnet)
library(survival)
source("moe/gate_features.R")

RD <- "moe/results"
RAW <- "moe/raw"
files <- list.files(RD, "rep_.*\\.rds$", full.names = TRUE)
cat(sprintf("Reading %d RDS files...\n", length(files)))

# Known feature names (fixed set, not dynamic)
FEATURES <- c("c_index","score_skew","score_kurt","score_var","score_q90_q10",
  "delta_c_index","prop_interact_sig","te_bin_var","te_slope","te_quadratic_p",
  "te_max_diff","te_bin_event_rate_range","censoring_rate","event_rate",
  "event_count_trt","event_count_ctrl","median_followup","mean_followup","n","p",
  "e_per_p","bootstrap_ci_sd","bootstrap_ci_mean","orig_pred_var","orig_pred_mean","orig_ambiguity")

# Config fields with prediction value (from generative model)
CONFIG_FEATURES <- c("n_predictive", "n_prognostic", "te_scale", "overlap", "b0",
  "prognostic_form", "corr", "te_start", "te_peak", "te_decay")

# Per-K prediction stats (computed by moe_simulation.R but not extracted)
K_STATS <- c("pred_var", "pred_mean", "ambiguity")

# New features computed from raw training data
NEW_FEATURES <- c("c_index_trt", "c_index_ctrl", "c_index_ratio",
  "te_int_max_z", "te_int_mean_z", "te_int_prop_sig",
  "corr_mean", "corr_max", "corr_prop_high")

rows <- list()
for (f in files) {
  r <- readRDS(f)
  row <- list(seed = r$seed,
    family = r$config$family,
    n_train = r$config$n_train,
    optimal_K = r$oracle_optimal_K,
    auc_K1 = r$aucs[1], auc_K2 = r$aucs[2], auc_K3 = r$aucs[3],
    auc_K4 = r$aucs[4], auc_K5 = r$aucs[5]
  )
  for (nm in FEATURES) {
    v <- r$features[[nm]]
    row[[nm]] <- if (is.null(v) || length(v) != 1) NA_real_ else v
  }
  # Config features
  for (nm in CONFIG_FEATURES) {
    v <- r$config[[nm]]
    row[[paste0("cfg_", nm)]] <- if (is.null(v) || length(v) != 1) NA_real_ else v
  }
  # Per-K prediction stats (K1..K5)
  for (K in 1:5) {
    for (stat in K_STATS) {
      key <- paste0("K", K, "_", stat)
      v <- r$features[[key]]
      row[[key]] <- if (is.null(v) || length(v) != 1) NA_real_ else v
    }
  }
  # New features from raw training data (within-arm, interactions, correlation)
  raw_path <- file.path(RAW, sprintf("%s_%d.rds", r$config$family, r$seed))
  if (file.exists(raw_path)) {
    d_raw <- readRDS(raw_path)
    train_raw <- d_raw$train
    if (!is.null(train_raw) && nrow(train_raw) > 20) {
      new_feats <- tryCatch(
        c(
          .extract_within_arm_features(train_raw),
          .extract_te_interaction_features(train_raw),
          .extract_correlation_features(train_raw)
        ),
        error = function(e) structure(rep(NA, 9), names = NEW_FEATURES)
      )
      for (nm in names(new_feats)) row[[nm]] <- new_feats[[nm]]
    }
  }
  rows[[length(rows) + 1]] <- row
  if (length(rows) %% 200 == 0) cat(sprintf("  %d/%d\n", length(rows), length(files)))
}

df <- rbindlist(rows, fill = TRUE)
fwrite(df, file.path(RD, "gate_training_data.csv"))
cat(sprintf("Written: %d rows x %d cols to moe/results/gate_training_data.csv\n", nrow(df), ncol(df)))
cat("Columns:", paste(names(df), collapse=", "), "\n")
