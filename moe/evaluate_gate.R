#!/usr/bin/env Rscript
# MoE-K Gate Evaluation: pre-trained gate vs real 5-fold CV
# Designed for the Omen (fast, sequential, Windows-compatible)
#
# Usage:
#   Rscript moe/evaluate_gate.R
#
# This script:
# 1. Reads gate training data (features + oracle optimal K)
# 2. Trains LASSO gate on 80% of data
# 3. For each test rep, runs real 5-fold CV over K=1..5
# 4. Compares: gate AUC vs CV AUC vs fixed K=4

library(glmnet)
library(data.table)
library(survival)

source("moe/scenario_generator.R")
source("R/stratified_cavboost.R")
source("R/rmst_cavboost_clean.R")

TAU <- 30
NR <- 30
.cov <- function(d) setdiff(names(d), c("trt01p","time","status","A","U","delta_tilde"))

auc_ <- function(p, l) {
  npos <- sum(l); nneg <- sum(!l)
  if (npos < 1 || nneg < 1) return(NA)
  r <- rank(p)
  (sum(r[as.logical(l)]) - npos * (npos + 1) / 2) / (npos * nneg)
}

# ---- Load gate training data ----
cat("Loading gate training data...\n")
d <- fread("moe/results/gate_training_data.csv")

feat_cols <- setdiff(names(d), c("family","n_train","optimal_K","auc_K1","auc_K2","auc_K3","auc_K4","auc_K5"))
cat(sprintf("  %d rows, %d features\n", nrow(d), length(feat_cols)))

X <- as.matrix(d[, ..feat_cols])
for (j in seq_len(ncol(X))) X[is.na(X[,j]), j] <- mean(X[,j], na.rm=TRUE)
X[!is.finite(X)] <- 0
y <- d$optimal_K

# ---- Train/test split ----
set.seed(20260530)
n <- nrow(X); idx <- sample(n, round(0.8 * n))
X_tr <- X[idx,]; y_tr <- y[idx]
X_te <- X[-idx,]; y_te <- y[-idx]
cat(sprintf("  Train: %d, Test: %d\n", length(y_tr), length(y_te)))

# ---- Train LASSO gate ----
cat("Training gate...\n")
cv_g <- cv.glmnet(X_tr, y_tr, family = "multinomial", alpha = 1, nfolds = 5)
gate <- glmnet(X_tr, y_tr, family = "multinomial", alpha = 1, lambda = cv_g$lambda.min)

# Gate predictions on test set
probs <- predict(gate, X_te, type = "response", s = cv_g$lambda.min)[,,1]
gate_K <- apply(probs, 1, which.max)

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
  set.seed(seed + 200)
  folds <- sample(rep(1:5, length.out = n_tr))
  cv_aucs <- numeric(5)

  for (K in 1:5) {
    fold_aucs <- numeric(5)
    for (f in 1:5) {
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
          set.seed(seed + K + f)
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
    cv_aucs[K] <- mean(fold_aucs[fold_aucs > 0], na.rm = TRUE)
  }

  # Get the actual test AUCs from the saved result
  test_aucs <- r$aucs
  oracle_opt <- r$oracle_optimal_K
  cv_K_opt <- which.max(cv_aucs)
  gate_K_val <- gate_K[i]

  results <- rbind(results, data.frame(
    family = cfg$family, n_train = cfg$n_train,
    oracle_auc = test_aucs[oracle_opt],
    gate_auc = if (!is.na(gate_K_val)) test_aucs[gate_K_val] else NA,
    cv_auc = if (!is.na(cv_K_opt)) test_aucs[cv_K_opt] else NA,
    k4_auc = test_aucs[4],
    oracle_K = oracle_opt,
    gate_K = gate_K_val,
    cv_K = cv_K_opt
  ))

  if (i %% 20 == 0) cat(sprintf("  %d/%d test reps done\n", i, length(test_files)))
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
cat(sprintf("%-18s %8.4f %8.4f %8.4f\n", "Fixed K=4", mean(results$k4_auc, na.rm = TRUE),
    mean(results$oracle_auc - results$k4_auc, na.rm = TRUE),
    max(results$oracle_auc - results$k4_auc, na.rm = TRUE)))

cat(sprintf("\nGate exact match: %.1f%%\n", 100 * mean(results$gate_K == results$oracle_K, na.rm = TRUE)))
cat(sprintf("CV exact match: %.1f%%\n", 100 * mean(results$cv_K == results$oracle_K, na.rm = TRUE)))

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
}

write.csv(results, "moe/results/gate_evaluation.csv", row.names = FALSE)
cat("\nResults saved to: moe/results/gate_evaluation.csv\n")