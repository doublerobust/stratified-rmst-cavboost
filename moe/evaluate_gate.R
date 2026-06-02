#!/usr/bin/env Rscript
# MoE-K Gate Evaluation: pre-trained gate vs real 5-fold CV
# Designed for the Omen (fast, sequential, Windows-compatible)
#
# Usage:
#   Rscript moe/evaluate_gate.R               # sequential
#   Rscript moe/evaluate_gate.R 0 10           # chunk 0 of 10 (parallel)
#   python3 moe/run_parallel_gate.py           # auto: 10 cores
#
# This script:
# 1. Reads gate training data (features + oracle optimal K)
# 2. Trains ranger RF gate on 80% of data
# 3. For each test rep, runs real 5-fold CV over K=1..5
# 4. Compares: gate AUC vs CV AUC vs fixed K=4

library(glmnet)
library(ranger)
library(data.table)
library(survival)

source("moe/scenario_generator.R")
source("R/stratified_cavboost.R")
source("R/rmst_cavboost_clean.R")
source("moe/gate_features.R")

TAU <- 30
NR <- 30
.cov <- function(d) setdiff(names(d), c("trt01p","time","status","A","U","delta_tilde"))

auc_ <- function(p, l) {
  npos <- sum(l); nneg <- sum(!l)
  if (npos < 1 || nneg < 1) return(NA)
  r <- rank(p)
  (sum(r[as.logical(l)]) - npos * (npos + 1) / 2) / (npos * nneg)
}

# ---- Parse parallel chunk args ----
args <- commandArgs(trailingOnly = TRUE)
chunk_idx <- 0L
n_chunks <- 1L
if (length(args) >= 2) {
  chunk_idx <- as.integer(args[1])
  n_chunks <- as.integer(args[2])
  cat(sprintf("Parallel mode: chunk %d of %d\n", chunk_idx + 1L, n_chunks))
}

# ---- Load gate training data ----
cat("Loading gate training data...\n")
d <- fread("moe/results/gate_training_data.csv")

exclude_cols <- c("seed","family","n_train","optimal_K","optimal_method","oracle_rate",
  "auc_K1","auc_K2","auc_K3","auc_K4","auc_K5","auc_VT",
  "cfg_n_predictive","cfg_n_prognostic","cfg_te_scale","cfg_b0",
  "cfg_prognostic_form","cfg_te_start","cfg_te_peak","cfg_te_decay",
  grep("^K[1-5]_", names(d), value = TRUE))
# Also exclude any non-numeric columns that slipped through
feat_cols <- setdiff(names(d), exclude_cols)
is_num <- sapply(d[, ..feat_cols], is.numeric)
feat_cols <- feat_cols[is_num]
cat(sprintf("  %d rows, %d features\n", nrow(d), length(feat_cols)))

X <- as.matrix(d[, ..feat_cols])
for (j in seq_len(ncol(X))) X[is.na(X[,j]), j] <- mean(X[,j], na.rm=TRUE)
X[!is.finite(X)] <- 0
y <- d$optimal_method

# ---- Train/test split ----
# Split at CONFIG level (not row level) to prevent leakage from
# per-rep data where multiple reps of the same config share a label
set.seed(20260530)
if (length(unique(paste0(d$family, d$n_train))) < nrow(d)) {
  # Per-rep data: group rows by config before splitting
  config_id <- rleidv(d[, c("family", "n_train")])
  unique_configs <- unique(config_id)
  n_configs <- length(unique_configs)
  train_configs <- sample(unique_configs, round(0.8 * n_configs))
  idx <- which(config_id %in% train_configs)
} else {
  # One row per config: simple row-level split is fine
  idx <- sample(nrow(X), round(0.8 * nrow(X)))
}
X_tr <- X[idx,]; y_tr <- y[idx]
X_te <- X[-idx,]; y_te <- y[-idx]
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

# ---- Real 5-fold CV on test set ----
cat(sprintf("Running real 5-fold CV on %d test reps...\n", length(y_te)))
cat("(This will take a while)\n\n")

files_all <- list.files("moe/results", "rep_.*\\.rds$", full.names = TRUE)
if (length(files_all) == 0) {
  cat("ERROR: No RDS files found in moe/results/\n")
  quit(status = 1)
}

d_test <- d[-idx, ]
test_files <- file.path("moe/results",
  paste0("rep_", d_test$family, "_", d_test$seed, ".rds"))
test_files <- test_files[file.exists(test_files)]

cat(sprintf("  Test files: %d\n", length(test_files)))

# ---- Parallel chunking ----
if (n_chunks > 1L) {
  chunk_size <- ceiling(length(test_files) / n_chunks)
  chunk_start <- chunk_idx * chunk_size + 1L
  chunk_end <- min((chunk_idx + 1L) * chunk_size, length(test_files))
  idx_chunk <- chunk_start:chunk_end
  test_files <- test_files[idx_chunk]
  d_test <- d_test[idx_chunk, , drop = FALSE]
  gate_K <- gate_K[idx_chunk]
  cat(sprintf("Chunk %d: files %d to %d (%d total)\n",
      chunk_idx + 1L, chunk_start, chunk_end, length(test_files)))
}

if (length(test_files) == 0) {
  cat("ERROR: gate_training_data.csv has no 'seed' column — re-run extract_gate_data.R\n")
  quit(status = 1)
}

results <- data.frame()
for (i in seq_along(test_files)) {
  r <- readRDS(test_files[i])
  cfg <- r$config; seed <- r$seed

  # Regenerate data (true oracle labels for CV)
  d_gen <- generate_scenario(family = cfg$family,
    n_predictive = cfg$n_predictive, n_prognostic = cfg$n_prognostic,
    te_scale = cfg$te_scale,
    overlap = cfg$overlap, b0 = cfg$b0,
    prognostic_form = cfg$prognostic_form,
    censoring_rate = cfg$censoring_rate,
    corr = cfg$corr, n_train = cfg$n_train, n_test = 2000,
    tau = TAU, seed = seed + 100)
  if (is.null(d_gen)) next

  train <- d_gen$train
  oracle <- as.logical(d_gen$train_oracle_label)
  if (all(oracle) || !any(oracle)) next  # degenerate fold
  zcols <- .cov(train)
  n_tr <- nrow(train)

  # 5-fold CV
  # Adaptive CV: try 5-fold -> 3-fold -> 2-fold, up to 3 random splits each
  run_cv_k <- function(train, oracle, zcols, K, base_seed) {
    for (n_folds in c(5L, 3L, 2L)) {
      for (attempt in 1:3) {
        seed_try <- base_seed + n_folds * 100 + attempt
        set.seed(seed_try)
        folds <- sample(rep(1:n_folds, length.out = nrow(train)))
        fold_aucs <- numeric(n_folds)
        for (f in 1:n_folds) {
          tr <- train[folds != f,]
          val <- train[folds == f,]
          oracle_val <- oracle[folds == f]
          if (sum(oracle_val) < 2 || sum(!oracle_val) < 2) next
          if (K == 1L) {
            fit <- tryCatch(train_rmst_cavboost(tr, tr$time, tr$status, TAU,
                              eta = 0.05, max_depth = 3, nr = NR),
                            error = function(e) NULL)
          } else {
            x_fold <- data.matrix(tr[, zcols, drop = FALSE])
            y_fold <- Surv(tr$time, tr$status)
            qq <- tryCatch({
              set.seed(seed_try + K + f)
              folds_inner <- sample(rep(1:5, length.out = nrow(tr)))
              lp <- numeric(nrow(tr))
              for (fi in 1:5) {
                tri <- which(folds_inner != fi); tei <- which(folds_inner == fi)
                cv_glm <- cv.glmnet(x_fold[tri,], y_fold[tri], family = "cox",
                                    alpha = 0.5, nfolds = 5, cox.ties = "breslow")
                lp[tei] <- drop(predict(cv_glm, x_fold[tei,], s = "lambda.min"))
              }
              unique(quantile(lp, seq(0, 1, 1/K), na.rm = TRUE))
            }, error = function(e) NULL)
            if (is.null(qq) || length(qq) < 2) next
            strata <- as.numeric(cut(lp, qq, include.lowest = TRUE, right = TRUE))
            fit <- tryCatch(train_stratified_cavboost(tr, tr$time, tr$status, TAU,
                                stratum = strata, eta = 0.05, max_depth = 3, nr = NR),
                            error = function(e) NULL)
          }
          if (is.null(fit)) next
          pred <- tryCatch(pred_subgroup(fit, val), error = function(e) NULL)
          if (!is.null(pred)) fold_aucs[f] <- auc_(pred, oracle_val)
        }
        fold_mean <- mean(fold_aucs[fold_aucs > 0], na.rm = TRUE)
        if (is.finite(fold_mean)) return(fold_mean)
      }
    }
    NA_real_
  }

  cv_aucs <- vapply(1:5, function(K) run_cv_k(train, oracle, zcols, K, seed + 200), numeric(1))
  
  # ---- CV for VT (ranger per arm, 5-fold) ----
  cv_vt_auc <- NA_real_
  zcols_vt <- setdiff(names(train), c("trt01p", "time", "status"))
  for (attempt in 1:3) {
    seed_try <- seed + 900 + attempt
    set.seed(seed_try)
    folds <- sample(rep(1:5, length.out = n_tr))
    fold_aucs <- numeric(5)
    for (f in 1:5) {
      tr <- train[folds != f,]
      val <- train[folds == f,]
      oracle_val <- oracle[folds == f]
      if (sum(oracle_val) < 2 || sum(!oracle_val) < 2) next
      ctrl <- tr[tr$trt01p == 0,]
      trt_d <- tr[tr$trt01p == 1,]
      if (nrow(ctrl) < 10 || nrow(trt_d) < 10) next
      rf_c <- tryCatch(ranger(Surv(time, status) ~ ., data = ctrl[, c("time", "status", zcols_vt)],
                           num.trees = 200, min.node.size = 10, seed = seed_try + f),
                     error = function(e) NULL)
      rf_t <- tryCatch(ranger(Surv(time, status) ~ ., data = trt_d[, c("time", "status", zcols_vt)],
                           num.trees = 200, min.node.size = 10, seed = seed_try + f + 100),
                     error = function(e) NULL)
      if (is.null(rf_c) || is.null(rf_t)) next
      pc <- predict(rf_c, val[, zcols_vt])
      pt <- predict(rf_t, val[, zcols_vt])
      tg <- seq(0, TAU, length.out = 200)
      surv_grid <- function(surv_mat, times, grid) {
        apply(surv_mat, 1, function(s) stats::approx(times, s, grid, rule = 2, yleft = 1)$y)
      }
      sc <- surv_grid(pc$survival, pc$unique.death.times, tg)
      st <- surv_grid(pt$survival, pt$unique.death.times, tg)
      vt_preds <- colMeans(st) * TAU - colMeans(sc) * TAU
      fold_aucs[f] <- auc_(vt_preds, oracle_val)
    }
    fold_mean <- mean(fold_aucs[fold_aucs > 0], na.rm = TRUE)
    if (is.finite(fold_mean)) { cv_vt_auc <- fold_mean; break }
  }

  # CV selects best among K1..K5 and VT
  cv_methods <- c(cv_aucs, cv_vt_auc)
  cv_method_opt <- if (length(cv_methods) > 0 && any(!is.na(cv_methods))) which.max(cv_methods) else NA_integer_

  # Get the actual test AUCs from the saved result
  test_aucs <- r$aucs
  oracle_opt <- if (!is.null(r$oracle_optimal_method)) r$oracle_optimal_method
                else if (!is.null(r$oracle_optimal_K)) r$oracle_optimal_K else NA_integer_
  cv_methods <- c(cv_aucs, cv_vt_auc)
  cv_method_opt <- if (length(cv_methods) > 0 && any(!is.na(cv_methods))) which.max(cv_methods) else NA_integer_
  gate_K_val <- gate_K[i]

  results <- rbind(results, data.frame(
    family = cfg$family, n_train = cfg$n_train,
    oracle_auc = test_aucs[oracle_opt],
    gate_auc = if (!is.na(gate_K_val)) test_aucs[gate_K_val] else NA,
    cv_auc = if (!is.na(cv_method_opt)) test_aucs[cv_method_opt] else NA,
    k4_auc = test_aucs[4],
    vt_auc = if (length(test_aucs) >= 6) test_aucs[6] else NA,
    oracle_method = oracle_opt,
    gate_method = gate_K_val,
    cv_method = cv_method_opt
  ))

  if (i %% 20 == 0) cat(sprintf("  %d/%d test reps done\n", i, length(test_files)))
  if (i %% 50 == 0) {
    dir.create("rds", showWarnings = FALSE)
    saveRDS(results, "rds/gate_intermediate_results.rds")
  }
}

# ---- Summary ----
cat("\n\n========== FINAL COMPARISON ==========\n\n")
cat(sprintf("%-18s %8s %8s %8s\n", "Method", "MeanAUC", "MeanGap", "MaxGap"))
cat(strrep("-", 42), "\n")
cat(sprintf("%-18s %8.4f %8s %8s\n", "Oracle", mean(results$oracle_auc, na.rm = TRUE), "-", "-"))
cat(sprintf("%-18s %8.4f %8.4f %8.4f\n", "Gate", mean(results$gate_auc, na.rm = TRUE),
    mean(results$oracle_auc - results$gate_auc, na.rm = TRUE),
    max(results$oracle_auc - results$gate_auc, na.rm = TRUE)))
cat(sprintf("%-18s %8.4f %8.4f %8.4f\n", "CV (real 5-fold)", mean(results$cv_auc, na.rm = TRUE),
    mean(results$oracle_auc - results$cv_auc, na.rm = TRUE),
    max(results$oracle_auc - results$cv_auc, na.rm = TRUE)))
cat(sprintf("%-18s %8.4f %8.4f %8.4f\n", "VT (always)", mean(results$vt_auc, na.rm = TRUE),
    mean(results$oracle_auc - results$vt_auc, na.rm = TRUE),
    max(results$oracle_auc - results$vt_auc, na.rm = TRUE)))
cat(sprintf("%-18s %8.4f %8.4f %8.4f\n", "Fixed K=4", mean(results$k4_auc, na.rm = TRUE),
    mean(results$oracle_auc - results$k4_auc, na.rm = TRUE),
    max(results$oracle_auc - results$k4_auc, na.rm = TRUE)))

cat(sprintf("\nGate exact match: %.1f%%\n", 100 * mean(results$gate_method == results$oracle_method, na.rm = TRUE)))
cat(sprintf("CV exact match: %.1f%%\n", 100 * mean(results$cv_method == results$oracle_method, na.rm = TRUE)))
cat(sprintf("Gate vs VT: Gate beats VT in %.1f%% of reps\n",
    100 * mean(results$gate_auc > results$vt_auc, na.rm = TRUE)))
cat(sprintf("VT vs CV:  VT beats CV in %.1f%% of reps\n",
    100 * mean(results$vt_auc > results$cv_auc, na.rm = TRUE)))
cat(sprintf("CV exact match: %.1f%%\n", 100 * mean(results$cv_method == results$oracle_method, na.rm = TRUE)))

cat("\n=== Method: VT-wins scenario ===\n")
vt_wins <- results[which(results$vt_auc > apply(results[,c("k4_auc","cv_auc","gate_auc")], 1, max, na.rm=TRUE)),]
if (nrow(vt_wins) > 0) {
  for (nn in c(200,300,500,1000)) {
    s <- vt_wins[vt_wins$n_train == nn,]
    if (nrow(s) < 2) next
    cat(sprintf("  n=%d (%d reps): VT wins when gate beats K=4 by %.4f\n",
        nn, nrow(s), mean(s$gate_auc - s$k4_auc, na.rm=TRUE)))
  }
}

cat("\n=== By n_train ===\n")
for (nn in c(200, 300, 500, 1000)) {
  s <- results[results$n_train == nn,]
  if (nrow(s) < 2) next
  cat(sprintf("\nn=%d (%d reps):\n", nn, nrow(s)))
  cat(sprintf("  Gate: %.4f (gap %.4f)\n", mean(s$gate_auc, na.rm = TRUE),
      mean(s$oracle_auc - s$gate_auc, na.rm = TRUE)))
  cat(sprintf("  CV:   %.4f (gap %.4f)\n", mean(s$cv_auc, na.rm = TRUE),
      mean(s$oracle_auc - s$cv_auc, na.rm = TRUE)))
  cat(sprintf("  K=4:  %.4f (gap %.4f)\n", mean(s$k4_auc, na.rm = TRUE),
      mean(s$oracle_auc - s$k4_auc, na.rm = TRUE)))
  cat(sprintf("  VT:   %.4f (gap %.4f)\n", mean(s$vt_auc, na.rm = TRUE),
      mean(s$oracle_auc - s$vt_auc, na.rm = TRUE)))
}

out_file <- if (n_chunks > 1L) {
  sprintf("moe/results/gate_evaluation_chunk_%d.csv", chunk_idx)
} else {
  "moe/results/gate_evaluation.csv"
}
write.csv(results, out_file, row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", out_file))