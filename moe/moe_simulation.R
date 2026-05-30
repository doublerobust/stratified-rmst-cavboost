#!/usr/bin/env Rscript
#'
#' MoE-K Simulation Pipeline
#' ==========================
#' Generates diverse scenarios, fits RMSTBoost at K=1..5,
#' extracts gate features, records oracle AUCs.
#'

suppressPackageStartupMessages({
  library(parallel)
})

# ---- Config ----
N_CONFIGS <- 1000
N_REPS <- 5
N_TOTAL <- N_CONFIGS * N_REPS
TAU <- 30
SEED_BASE <- 20260601
NR <- 30  # boosting iterations (reduced from 50 — sufficient for K ranking)

# Use the directory containing this script as the base
REPO_DIR <- getwd()
MOE_DIR <- file.path(REPO_DIR, "moe")
RAW_DIR <- file.path(MOE_DIR, "raw")
RESULTS_DIR <- file.path(MOE_DIR, "results")

dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Source project files ----
source(file.path(MOE_DIR, "scenario_generator.R"))
source(file.path(MOE_DIR, "gate_features.R"))
source(file.path(REPO_DIR, "R", "stratified_cavboost.R"))
source(file.path(REPO_DIR, "R", "rmst_cavboost_clean.R"))

# ---- AUC helper ----
auc_ <- function(p, l) {
  npos <- sum(l); nneg <- sum(!l)
  if (npos < 1 || nneg < 1) return(NA)
  r <- rank(p)
  (sum(r[as.logical(l)]) - npos * (npos + 1) / 2) / (npos * nneg)
}

# ---- Covariate column names ----
.covariate_cols <- function(dat) {
  setdiff(names(dat), c("trt01p", "time", "status", "A", "U", "delta_tilde"))
}

# ---- Process one rep ----
process_rep <- function(config_idx, all_configs) {
  cfg <- all_configs[config_idx, ]

  cat(sprintf("[config %d] starting (family=%s, n_train=%d)\n",
              config_idx, cfg$family, cfg$n_train))
  utils::flush.console()

  for (rep in seq_len(N_REPS)) {
    seed <- cfg$seed + rep

    # Skip if already done
    result_path <- file.path(RESULTS_DIR, sprintf("rep_%s_%d.rds", cfg$family, seed))
    if (file.exists(result_path)) {
      cat(sprintf("[config %d] rep %d skip (exists)\n", config_idx, rep))
      utils::flush.console()
      next
    }

    # Generate scenario
    cat(sprintf("[config %d] rep %d generate_scenario...\n", config_idx, rep))
    utils::flush.console()
    d <- generate_scenario(
      family = cfg$family,
      n_predictive = cfg$n_predictive,
      n_prognostic = cfg$n_prognostic,
      overlap = cfg$overlap,
      b0 = cfg$b0,
      prognostic_form = cfg$prognostic_form,
      censoring_rate = cfg$censoring_rate,
      corr = cfg$corr,
      n_train = cfg$n_train,
      n_test = 2000,
      tau = TAU,
      seed = seed,
      save_dir = RAW_DIR
    )

    if (is.null(d)) {
      cat(sprintf("[config %d] rep %d invalid (null)\n", config_idx, rep))
      utils::flush.console()
      next
    }

    train <- d$train
    test <- d$test
    oracle <- d$oracle_label
    zcols <- .covariate_cols(train)

    # ---- Fit RMSTBoost for K = 1..5 ----
    k_values <- 1:5
    aucs <- setNames(rep(NA_real_, 5), paste0("auc_K", k_values))
    preds_list <- list()

    for (K in k_values) {
      cat(sprintf("[config %d] rep %d K=%d...\n", config_idx, rep, K))
      utils::flush.console()

      if (K == 1L) {
        fit <- tryCatch(
          train_rmst_cavboost(train, train$time, train$status, TAU,
                              eta = 0.05, max_depth = 3, nr = NR),
          error = function(e) NULL
        )
      } else {
        strata <- tryCatch(
          crossfit_prognostic_strata(train, zcols, K = K, seed = seed + K),
          error = function(e) NULL
        )
        if (is.null(strata) || length(unique(strata)) < 2) {
          cat(sprintf("[config %d] rep %d K=%d strata invalid\n", config_idx, rep, K))
          utils::flush.console()
          next
        }

        fit <- tryCatch(
          train_stratified_cavboost(train, train$time, train$status, TAU,
                                    stratum = strata,
                                    eta = 0.05, max_depth = 3, nr = NR),
          error = function(e) NULL
        )
      }

      if (!is.null(fit)) {
        cat(sprintf("[config %d] rep %d K=%d pred_subgroup...\n", config_idx, rep, K))
        utils::flush.console()
        pred <- tryCatch(pred_subgroup(fit, test), error = function(e) NULL)
        if (!is.null(pred)) {
          preds_list[[K]] <- pred
          aucs[K] <- auc_(pred, oracle)
          cat(sprintf("[config %d] rep %d K=%d AUC=%.3f\n", config_idx, rep, K, aucs[K]))
          utils::flush.console()
        }
      }
    }

    # ---- Compute gate features from training data ----
    # Fit base Orig model for internal behavior features
    cat(sprintf("[config %d] rep %d fit_orig...\n", config_idx, rep))
    utils::flush.console()
    fit_orig <- tryCatch(
      train_rmst_cavboost(train, train$time, train$status, TAU,
                          eta = 0.05, max_depth = 3, nr = NR),
      error = function(e) NULL
    )
    train_preds <- if (!is.null(fit_orig)) {
      tryCatch(pred_subgroup(fit_orig, train), error = function(e) NULL)
    } else NULL

    cat(sprintf("[config %d] rep %d compute_gate_features...\n", config_idx, rep))
    utils::flush.console()
    features <- tryCatch(
      compute_gate_features(train, tau = TAU, fit_orig = fit_orig, model_preds = train_preds),
      error = function(e) structure(rep(NA_real_, 30), names = paste0("feat_", 1:30))
    )

    # Add per-K prediction stats (only if we have predictions)
    if (length(preds_list) > 0) {
      for (K in k_values) {
        if (!is.null(preds_list[[K]])) {
          pk <- preds_list[[K]]
          features[paste0("K", K, "_pred_var")] <- var(pk, na.rm = TRUE)
          features[paste0("K", K, "_pred_mean")] <- mean(pk, na.rm = TRUE)
          features[paste0("K", K, "_ambiguity")] <- mean(pk > 0.4 & pk < 0.6, na.rm = TRUE)
        }
      }
    }

    # ---- Save result ----
    result <- list(
      config_idx = config_idx,
      config = cfg,
      rep = rep,
      seed = seed,
      aucs = aucs,
      oracle_optimal_K = if (length(preds_list) > 0) which.max(aucs) else NA_integer_,
      features = features,
      oracle_rate = mean(oracle)
    )

    cat(sprintf("[config %d] rep %d saveRDS...\n", config_idx, rep))
    utils::flush.console()
    saveRDS(result, result_path)
    cat(sprintf("[config %d] rep %d done\n", config_idx, rep))
    utils::flush.console()
  }

  cat(sprintf("[config %d] complete\n", config_idx))
  utils::flush.console()
  invisible(NULL)
}

# ---- Main ----
cat(sprintf("MoE-K Simulation Pipeline\n"))
cat(sprintf("Configs: %d, Reps per config: %d, Total: %d\n", N_CONFIGS, N_REPS, N_TOTAL))
cat(sprintf("Raw data: %s\n", RAW_DIR))
cat(sprintf("Results: %s\n", RESULTS_DIR))
cat("========================================\n\n")

set.seed(SEED_BASE)
configs <- sample_configurations(N_CONFIGS)
cat(sprintf("Generated %d scenario configurations\n", nrow(configs)))
cat(sprintf("Family distribution:\n"))
print(table(configs$family))
cat(sprintf("Train sample size distribution:\n"))
print(table(configs$n_train))
cat("\n")

cat(sprintf("Processing sequentially (Windows compatibility)...\n\n"))

results <- lapply(seq_len(N_CONFIGS), process_rep,
                  all_configs = configs)

cat("\n=== Simulation Complete ===\n")

# ---- Summarize ----
result_files <- list.files(RESULTS_DIR, "rep_.*\\.rds$", full.names = TRUE)
cat(sprintf("Files written: %d / %d\n", length(result_files), N_TOTAL))

if (length(result_files) > 0) {
  summary_list <- lapply(result_files, function(f) {
    r <- readRDS(f)
    data.frame(
      family = r$config$family,
      n_train = r$config$n_train,
      auc_K1 = r$aucs[1], auc_K2 = r$aucs[2], auc_K3 = r$aucs[3],
      auc_K4 = r$aucs[4], auc_K5 = r$aucs[5],
      oracle_optimal_K = r$oracle_optimal_K,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_list)

  cat("\nOptimal K distribution (overall):\n")
  print(table(summary_df$oracle_optimal_K))

  cat("\nOptimal K by family:\n")
  tbl <- table(summary_df$family, summary_df$oracle_optimal_K)
  tbl_prop <- prop.table(tbl, 1)
  print(round(tbl_prop * 100, 1))

  cat("\nOptimal K by n_train:\n")
  tbl_n <- table(summary_df$n_train, summary_df$oracle_optimal_K)
  tbl_n_prop <- prop.table(tbl_n, 1)
  print(round(tbl_n_prop * 100, 1))

  # AUC improvement vs K=4
  summary_df$auc_improvement <- summary_df$auc_K4 - apply(summary_df[, paste0("auc_K", 1:5)], 1, max, na.rm = TRUE)
  cat(sprintf("\nAUC loss from using K=4 vs optimal: mean = %.4f, median = %.4f\n",
              mean(abs(summary_df$auc_improvement), na.rm = TRUE),
              median(abs(summary_df$auc_improvement), na.rm = TRUE)))

  write.csv(summary_df, file.path(RESULTS_DIR, "summary.csv"), row.names = FALSE)
  cat("Summary saved to:", file.path(RESULTS_DIR, "summary.csv"), "\n")
}

cat("\nDone.\n")