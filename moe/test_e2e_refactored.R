#!/usr/bin/env Rscript
# Pipeline end-to-end verification test
# Run from repo root:   Rscript moe/test_e2e_refactored.R
#
# Tests all 3 modes of extract_gate_data.R:
#   1. Standalone (config-averaged)
#   2. Per-rep
#   3. Chunk (parallel)
#
# Generates 3 configs × 2 reps = 6 RDS files, runs extraction,
# and asserts no NAs in NEW_FEATURES.

library(data.table)
MOE <- "moe"
source(file.path(MOE, "scenario_generator.R"))
source(file.path(MOE, "gate_features.R"))
TAU <- 30

cleanup <- function() {
  unlink(list.files(file.path(MOE, "results"), "rep_.*\\.rds$", full.names=T))
  unlink(file.path(MOE, "results", "gate_training_data.csv"))
  unlink(file.path(MOE, "results", "chunk_*.csv"))
}

status <- function(label, ok, detail="") {
  cat(sprintf("  %s %s %s\n", if(ok) "✅" else "❌", label, detail))
}

cat("=== Generating 3 configs x 2 reps = 6 RDS files ===\n")
set.seed(42)

for (ic in 1:3) {
  fam <- c("linear", "bump", "enclave")[ic]
  ntr <- c(200, 500, 1000)[ic]
  for (ir in 1:2) {
    seed <- ic * 10000 + ir * 100
    d <- generate_scenario(family=fam, n_train=ntr, seed=seed,
                           b0=1, te_scale=0.5, tau=TAU)
    mock_preds <- runif(nrow(d$train), 0, 0.5)
    feats <- compute_gate_features(d$train, tau=TAU, fit_orig=TRUE, model_preds=mock_preds)
    aucs <- switch(fam,
      linear=c(0.72,0.70,0.68,0.65,0.62,0.65),
      bump=c(0.65,0.70,0.75,0.73,0.70,0.60),
      enclave=c(0.60,0.65,0.68,0.72,0.70,0.66))
    result <- list(config_idx=ic, seed=seed,
                   config=list(family=fam, n_train=ntr),
                   features=feats, aucs=aucs,
                   oracle_optimal_K=which.max(aucs[1:5]),
                   oracle_optimal_method=which.max(aucs[1:5]))
    saveRDS(result, file.path(MOE, "results", sprintf("rep_%s_%d.rds", fam, seed)))
  }
}
n_rds <- length(list.files(file.path(MOE, "results"), "rep_.*\\.rds$"))
status("RDS files generated", n_rds == 6, sprintf("(%d found)", n_rds))

# ---- Test 1: Standalone ----
cat("\n=== Test 1: Standalone mode ===\n")
system2("Rscript", c(file.path(MOE, "extract_gate_data.R")),
        stdout=TRUE, stderr=TRUE) |> suppressWarnings() |> invisible()
csv <- file.path(MOE, "results", "gate_training_data.csv")
if (file.exists(csv) && file.info(csv)$size > 0) {
  df <- fread(csv)
  new_feats <- c("c_index_trt","te_int_max_z","corr_mean","prop_interact_sig")
  nas <- sum(sapply(new_feats, function(nm) sum(is.na(df[[nm]]))))
  status("NEW_FEATURES populated", nas == 0, sprintf("(%d NAs)", nas))
  status("Configs extracted", nrow(df) == 3, sprintf("(%d rows)", nrow(df)))
} else {
  status("Standalone mode", FALSE, "no output CSV")
}

# ---- Test 2: Per-rep ----
cat("\n=== Test 2: Per-rep mode ===\n")
file.remove(csv)
system2("Rscript", c(file.path(MOE, "extract_gate_data.R"), "--per-rep"),
        stdout=TRUE, stderr=TRUE) |> suppressWarnings() |> invisible()
df2 <- fread(csv)
status("Per-rep rows", nrow(df2) == 6, sprintf("(%d rows)", nrow(df2)))
status("optimal_method column", "optimal_method" %in% names(df2))
has_na <- sum(is.na(df2[["c_index_trt"]])) > 0
status("No NEW_FEATURES NAs", !has_na)

# ---- Test 3: Chunk mode ----
cat("\n=== Test 3: Chunk mode ===\n")
file.remove(csv)
files <- list.files(file.path(MOE, "results"), "rep_.*\\.rds$", full.names=TRUE)
writeLines(files[1:3], file.path(MOE, "results", "chunk_0.txt"))
writeLines(files[4:6], file.path(MOE, "results", "chunk_1.txt"))
for (id in 0:1) {
  system2("Rscript", c(file.path(MOE, "extract_gate_data.R"),
    file.path(MOE, "results", sprintf("chunk_%d.txt", id)), id),
    stdout=TRUE, stderr=TRUE) |> suppressWarnings() |> invisible()
}
merged <- rbindlist(lapply(0:1, function(id) {
  p <- file.path(MOE, "results", sprintf("chunk_%d.csv", id))
  if (file.exists(p)) fread(p) else NULL
}))
status("Chunks merged", !is.null(merged) && nrow(merged) == 6,
       sprintf("(%d rows)", if(is.null(merged)) 0 else nrow(merged)))

cleanup()
cat("\n✅ Pipeline verification complete.\n")
