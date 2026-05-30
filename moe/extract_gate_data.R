#!/usr/bin/env Rscript
# Extract gate training data from RDS results
# Dumps features + optimal K to CSV (since RDS is gitignored)
library(data.table)
library(glmnet)

RD <- "moe/results"
files <- list.files(RD, "rep_.*\\.rds$", full.names = TRUE)
cat(sprintf("Reading %d RDS files...\n", length(files)))

# Known feature names (fixed set, not dynamic)
FEATURES <- c("c_index","score_skew","score_kurt","score_var","score_q90_q10",
  "delta_c_index","prop_interact_sig","te_bin_var","te_slope","te_quadratic_p",
  "te_max_diff","te_bin_event_rate_range","censoring_rate","event_rate",
  "event_count_trt","event_count_ctrl","median_followup","mean_followup","n","p",
  "e_per_p","bootstrap_ci_sd","bootstrap_ci_mean","orig_pred_var","orig_pred_mean","orig_ambiguity")

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
  rows[[length(rows) + 1]] <- row
}

df <- rbindlist(rows, fill = TRUE)
fwrite(df, file.path(RD, "gate_training_data.csv"))
cat(sprintf("Written: %d rows x %d cols to moe/results/gate_training_data.csv\n", nrow(df), ncol(df)))
cat("Columns:", paste(names(df), collapse=", "), "\n")
