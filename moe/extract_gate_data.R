#!/usr/bin/env Rscript
# Extract gate training data from RDS results + raw data
# Dumps features + optimal K + richer meta-features to CSV
#
# Standalone:
#   Rscript moe/extract_gate_data.R                          # all files, averaged by config
#   Rscript moe/extract_gate_data.R --per-rep                # all files, keep per-rep rows
#
# Chunk mode (for parallel extraction):
#   Rscript moe/extract_gate_data.R chunk_file chunk_id      # process one chunk, write chunk_N.csv

library(data.table)
library(glmnet)
library(survival)
source("moe/gate_features.R")

RD <- "moe/results"
RAW <- "moe/raw"
args <- commandArgs(TRUE)

# ── Detect mode ─────────────────────────────────────────
CHUNK_MODE <- length(args) >= 2 && file.exists(args[1])
PER_REP <- "--per-rep" %in% args

if (CHUNK_MODE) {
  # Chunk mode: args = [chunk_file, chunk_id]
  chunk_file <- args[1]
  chunk_id <- as.integer(args[2])
  files <- readLines(chunk_file)
  out_csv <- file.path(RD, sprintf("chunk_%d.csv", chunk_id))
} else {
  files <- list.files(RD, "rep_.*\\.rds$", full.names = TRUE)
  out_csv <- file.path(RD, "gate_training_data.csv")
}

cat(sprintf("Reading %d RDS files...\n", length(files)))

# ── Feature lists ───────────────────────────────────────
FEATURES <- c("c_index","score_skew","score_kurt","score_var","score_q90_q10",
  "delta_c_index","te_bin_var","te_slope","te_quadratic_p",
  "te_max_diff","te_bin_event_rate_range","censoring_rate","event_rate",
  "event_count_trt","event_count_ctrl","median_followup","mean_followup","n","p",
  "e_per_p","bootstrap_ci_sd","bootstrap_ci_mean","orig_pred_var","orig_pred_mean","orig_ambiguity")

CONFIG_FEATURES <- c("n_predictive", "n_prognostic", "te_scale", "overlap", "b0",
  "prognostic_form", "corr", "te_start", "te_peak", "te_decay")

NEW_FEATURES <- c("c_index_trt", "c_index_ctrl", "c_index_ratio",
  "te_int_max_z", "te_int_mean_z", "te_int_prop_sig",
  "corr_mean", "corr_max", "corr_prop_high",
  "prop_interact_sig", "trt_main_p")

# ── Extract one file ────────────────────────────────────
extract_one <- function(f) {
  r <- readRDS(f)
  row <- list(
    config_idx = r$config_idx,
    seed = r$seed,
    family = r$config$family,
    n_train = r$config$n_train,
    auc_K1 = r$aucs[1], auc_K2 = r$aucs[2], auc_K3 = r$aucs[3],
    auc_K4 = r$aucs[4], auc_K5 = r$aucs[5],
    auc_VT = if (length(r$aucs) >= 6) r$aucs[6] else NA
  )
  
  # FEATURES from RDS
  for (nm in FEATURES) {
    if (!nm %in% names(r$features)) { row[[nm]] <- NA_real_; next }
    v <- r$features[[nm]]
    row[[nm]] <- if (is.null(v) || length(v) != 1) NA_real_ else v
  }
  
  # Config features
  for (nm in CONFIG_FEATURES) {
    v <- r$config[[nm]]
    row[[paste0("cfg_", nm)]] <- if (is.null(v) || length(v) != 1) NA_real_ else v
  }
  
  # NEW_FEATURES from raw data
  raw_path <- file.path(RAW, sprintf("%s_%d.rds", r$config$family, r$seed))
  if (file.exists(raw_path)) {
    d_raw <- readRDS(raw_path)
    train_raw <- d_raw$train
    if (!is.null(train_raw) && nrow(train_raw) > 20) {
      new_feats <- tryCatch(c(
        .extract_within_arm_features(train_raw),
        .extract_interaction_features(train_raw),
        .extract_te_interaction_features(train_raw),
        .extract_correlation_features(train_raw)
      ), error = function(e)
        structure(rep(NA, length(NEW_FEATURES)), names = NEW_FEATURES))
      for (nm in names(new_feats)) row[[nm]] <- new_feats[[nm]]
    }
  }
  
  row
}

# ── Process all files ───────────────────────────────────
rows <- list()
for (i in seq_along(files)) {
  rows[[i]] <- extract_one(files[i])
  if (i %% 200 == 0) cat(sprintf("  %d/%d\n", i, length(files)))
}

df <- rbindlist(rows, fill = TRUE)

# ── Post-processing ────────────────────────────────────
if (!CHUNK_MODE && "config_idx" %in% names(df)) {
  cat(sprintf("Aggregating %d rows by config_idx...\n", nrow(df)))
  
  auc_cols <- c("auc_K1", "auc_K2", "auc_K3", "auc_K4", "auc_K5", "auc_VT")
  exclude_from_avg <- c("config_idx", "seed", "family", "n_train",
    "auc_cols", "cfg_prognostic_form", "cfg_te_start", "cfg_te_peak", "cfg_te_decay")
  feat_cols <- setdiff(names(df), exclude_from_avg)
  is_num <- sapply(df[, ..feat_cols], is.numeric)
  feat_cols <- feat_cols[is_num]
  
  # Average numeric columns per config
  agg <- df[, lapply(.SD, mean, na.rm = TRUE),
            by = config_idx, .SDcols = c(auc_cols, feat_cols)]
  
  # Add back non-numeric config-level columns
  first_cols <- intersect(c("family", "n_train", "cfg_prognostic_form"), names(df))
  first_vals <- unique(df[, c("config_idx", first_cols), with = FALSE])
  agg <- merge(agg, first_vals, by = "config_idx")
  
  # Oracle: best K by mean AUC
  auc_mat <- as.matrix(agg[, ..auc_cols])
  na_all <- rowSums(is.na(auc_mat)) == ncol(auc_mat)
  if (any(na_all)) agg <- agg[!na_all, ]
  auc_mat <- as.matrix(agg[, ..auc_cols])
  auc_mat[is.na(auc_mat)] <- -Inf
  agg$optimal_method <- max.col(auc_mat, ties.method = "first")
  
  # Keep a reference seed
  first_seeds <- df[, .(seed = seed[1]), by = config_idx]
  agg <- merge(agg, first_seeds, by = "config_idx")
  agg[, config_idx := NULL]
  
  if (PER_REP) {
    # Keep per-rep rows but add config-level oracle label
    df <- merge(df, agg[, .(config_idx = unique(df$config_idx),
                            optimal_method = agg$optimal_method),
                        by = config_idx], by = "config_idx")
    # Keep all per-rep columns
    df[, config_idx := NULL]
    fwrite(df, out_csv)
  } else {
    df <- agg
    fwrite(df, out_csv)
  }
} else {
  fwrite(df, out_csv)
}

cat(sprintf("Written: %d rows x %d cols to %s\n", nrow(df), ncol(df), out_csv))
cat("Done.\n")
