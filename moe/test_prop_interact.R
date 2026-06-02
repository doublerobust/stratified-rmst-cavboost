#!/usr/bin/env Rscript
# Quick test: verify prop_interact_sig is correctly extracted
source("moe/gate_features.R")
library(data.table)

# 1. Check if feature exists in simulation results
files <- list.files("moe/results", "rep_.*\\.rds$", full.names = TRUE)
cat(sprintf("Total result files: %d\n\n", length(files)))

# Check first 20 files
n_check <- min(20, length(files))
for (i in seq_len(n_check)) {
  r <- readRDS(files[i])
  val <- r$features[["prop_interact_sig"]]
  cat(sprintf("Result %d (family=%s): prop_interact_sig=%s\n",
              i, r$config$family, if (is.null(val)) "NULL" else if (is.na(val)) "NA" else round(val, 4)))
}

cat("\n--- Checking raw data recomputation ---\n")

# 2. Check if raw data files exist for first 5 results
for (i in seq_len(min(5, length(files)))) {
  r <- readRDS(files[i])
  raw_path <- file.path("moe/raw", sprintf("%s_%d.rds", r$config$family, r$seed))
  
  if (!file.exists(raw_path)) {
    cat(sprintf("Result %d: RAW file missing: %s\n", i, raw_path))
    next
  }
  
  d_raw <- readRDS(raw_path)
  train_raw <- d_raw$train
  cat(sprintf("Result %d (family=%s, n=%d): raw file found\n", i, r$config$family, nrow(train_raw)))
  
  # Compute prop_interact_sig from raw training data
  feats <- tryCatch(
    .extract_interaction_features(train_raw),
    error = function(e) c(prop_interact_sig = NA)
  )
  cat(sprintf("  prop_interact_sig from raw data: %s\n",
              if (is.na(feats["prop_interact_sig"])) "NA" else round(feats["prop_interact_sig"], 4)))
  cat(sprintf("  delta_c_index: %s\n",
              if (is.na(feats["delta_c_index"])) "NA" else round(feats["delta_c_index"], 4)))
  cat(sprintf("  trt_main_p: %s\n",
              if (is.na(feats["trt_main_p"])) "NA" else round(feats["trt_main_p"], 4)))
  cat("\n")
}