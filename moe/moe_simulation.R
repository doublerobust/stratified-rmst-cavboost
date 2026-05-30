#!/usr/bin/env Rscript
#'
#' MoE-K Simulation Pipeline
#' ==========================
#' Generates diverse scenarios, fits RMSTBoost at K=1..5,
#' extracts gate features, records oracle AUCs.
#' 

suppressPackageStartupMessages({
  library(parallel)
  library(survival)
  library(glmnet)
})

# ---- Config (defaults; override via N_CONFIGS <<- before source()) ----
if (!exists("N_CONFIGS")) N_CONFIGS <- 1000
if (!exists("N_REPS")) N_REPS <- 5
if (!exists("N_WORKERS")) N_WORKERS <- 11
N_TOTAL <- N_CONFIGS * N_REPS
TAU <- 30
SEED_BASE <- 20260601
NR <- 30  # boosting iterations
NTHREAD <- 2  # XGBoost threads per worker (N_WORKERS * NTHREAD <= total cores)
PARALLEL <- TRUE

# Use the directory containing this script as the base
REPO_DIR <- getwd()
MOE_DIR <- file.path(REPO_DIR, "moe")
RAW_DIR <- file.path(MOE_DIR, "raw")
RESULTS_DIR <- file.path(MOE_DIR, "results")

dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Source project files (top-level: needed for sample_configurations, etc.) ----
source(file.path(MOE_DIR, "scenario_generator.R"))
source(file.path(MOE_DIR, "gate_features.R"))
source(file.path(REPO_DIR, "R", "stratified_cavboost.R"))
source(file.path(REPO_DIR, "R", "rmst_cavboost_clean.R"))

# ---- Process one rep ----
process_rep <- function(config_idx, all_configs, moe_dir, repo_dir,
                         raw_dir, results_dir, tau, nr, nthread) {
  # Source project files inside each worker
  source(file.path(moe_dir, "scenario_generator.R"))
  source(file.path(moe_dir, "gate_features.R"))
  source(file.path(repo_dir, "R", "stratified_cavboost.R"))
  source(file.path(repo_dir, "R", "rmst_cavboost_clean.R"))

  # Helper functions (defined inside body so PSOCK workers inherit them)
  auc_ <- function(p, l) {
    npos <- sum(l); nneg <- sum(!l)
    if (npos < 1 || nneg < 1) return(NA)
    r <- rank(p)
    (sum(r[as.logical(l)]) - npos * (npos + 1) / 2) / (npos * nneg)
  }
  
  .covariate_cols <- function(dat) {
    setdiff(names(dat), c("trt01p", "time", "status", "A", "U", "delta_tilde"))
  }

  cfg <- all_configs[config_idx, ]

  for (rep in seq_len(N_REPS)) {
    seed <- cfg$seed + rep

    # Skip if already done
    result_path <- file.path(results_dir, sprintf("rep_%s_%d.rds", cfg$family, seed))
    if (file.exists(result_path)) next

    # Generate scenario
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
      tau = tau,
      seed = seed,
      save_dir = raw_dir
    )

    if (is.null(d)) next

    train <- d$train
    test <- d$test
    oracle <- d$oracle_label
    zcols <- .covariate_cols(train)

    # ---- Fit RMSTBoost for K = 1..5 ----
    k_values <- 1:5
    aucs <- setNames(rep(NA_real_, 5), paste0("auc_K", k_values))
    preds_list <- vector("list", 5)

    # ---- Compute prognostic score ONCE ----
    prog_lp <- NULL
    for (attempt in 1:5) {
      prog_lp <- tryCatch({
        n_tr <- nrow(train)
        x <- data.matrix(train[, zcols, drop = FALSE])
        y <- survival::Surv(train$time, train$status)
        set.seed(seed + 999 + attempt * 100)
        folds <- sample(rep(1:5, length.out = n_tr))
        lp <- numeric(n_tr)
        for (fold in 1:5) {
          tr_idx <- which(folds != fold); te_idx <- which(folds == fold)
          cv <- suppressWarnings(glmnet::cv.glmnet(
            x[tr_idx, , drop = FALSE], y[tr_idx],
            family = "cox", alpha = 0.5, nfolds = 5, 
            cox.ties = "breslow"))
          lp[te_idx] <- drop(stats::predict(cv, x[te_idx, , drop = FALSE], s = "lambda.min"))
        }
        lp
      }, error = function(e) NULL)
      if (!is.null(prog_lp)) break
    }

    for (K in k_values) {
      if (K == 1L) {
        fit <- tryCatch(
          train_rmst_cavboost(train, train$time, train$status, tau,
                              eta = 0.05, max_depth = 3, nr = nr,
                              nthread = nthread),
          error = function(e) NULL
        )
      } else if (!is.null(prog_lp)) {
        qq <- unique(stats::quantile(prog_lp, seq(0, 1, 1 / K), na.rm = TRUE))
        if (length(qq) >= 2) {
          strata <- as.numeric(cut(prog_lp, qq, include.lowest = TRUE, right = TRUE))
          fit <- tryCatch(
            train_stratified_cavboost(train, train$time, train$status, tau,
                                      stratum = strata,
                                      eta = 0.05, max_depth = 3, nr = nr,
                                      nthread = nthread),
            error = function(e) NULL
          )
        } else {
          fit <- NULL
        }
      } else {
        fit <- NULL
      }

      if (!is.null(fit)) {
        pred <- tryCatch(pred_subgroup(fit, test), error = function(e) NULL)
        if (!is.null(pred)) {
          preds_list[[K]] <- pred
          aucs[K] <- auc_(pred, oracle)
        }
      }
    }

    # ---- Compute gate features ----
    fit_orig <- tryCatch(
      train_rmst_cavboost(train, train$time, train$status, tau,
                          eta = 0.05, max_depth = 3, nr = nr,
                          nthread = nthread),
      error = function(e) NULL
    )
    train_preds <- if (!is.null(fit_orig)) {
      tryCatch(pred_subgroup(fit_orig, train), error = function(e) NULL)
    } else NULL

    features <- tryCatch(
      compute_gate_features(train, tau = tau, fit_orig = fit_orig, model_preds = train_preds),
      error = function(e) structure(rep(NA_real_, 30), names = paste0("feat_", 1:30))
    )

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

    saveRDS(result, result_path)
  }

  if (config_idx %% 20 == 0) {
    cat(sprintf("  progress: %d / %d configs\n", config_idx, N_CONFIGS))
    utils::flush.console()
  }

  gc()
  invisible(NULL)
}

# ---- Main ----
cat(sprintf("MoE-K Simulation Pipeline\n"))
cat(sprintf("Configs: %d, Reps per config: %d, Total: %d\n", N_CONFIGS, N_REPS, N_TOTAL))
cat(sprintf("Workers: %d, XGBoost threads per worker: %d\n", if (PARALLEL) N_WORKERS else 1, NTHREAD))
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

if (PARALLEL && N_CONFIGS >= 10) {
  cat(sprintf("Starting parallel processing with %d workers...\n\n", N_WORKERS))
  utils::flush.console()

  cl <- makePSOCKcluster(N_WORKERS, outfile = "")

  # Export variables needed by workers
  clusterExport(cl, c("N_REPS"), envir = environment())

  # Process in parallel
  parLapply(cl, seq_len(N_CONFIGS), process_rep,
            all_configs = configs,
            moe_dir = MOE_DIR,
            repo_dir = REPO_DIR,
            raw_dir = RAW_DIR,
            results_dir = RESULTS_DIR,
            tau = TAU,
            nr = NR,
            nthread = NTHREAD)

  stopCluster(cl)
} else {
  cat("Processing sequentially...\n\n")
  utils::flush.console()

  lapply(seq_len(N_CONFIGS), process_rep,
         all_configs = configs,
         moe_dir = MOE_DIR,
         repo_dir = REPO_DIR,
         raw_dir = RAW_DIR,
         results_dir = RESULTS_DIR,
         tau = TAU,
         nr = NR,
         nthread = 1L)
}

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

  # Flag near-perfect AUCs for diagnostic
  perfect <- apply(summary_df[, paste0("auc_K", 1:5)] > 0.999, 1, any, na.rm = TRUE)
  n_perfect <- sum(perfect, na.rm = TRUE)
  if (n_perfect > 0) {
    cat(sprintf("\n\u26a0\ufe0f  %d / %d reps have AUC > 0.999 (near-perfect prediction)\n",
                n_perfect, nrow(summary_df)))
    perfect_df <- summary_df[perfect, ]
    cat("\nBy family:\n")
    print(table(perfect_df$family))
    cat("\nBy n_train:\n")
    print(table(perfect_df$n_train))
  } else {
    cat("\nNo near-perfect AUCs detected.\n")
  }

  write.csv(summary_df, file.path(RESULTS_DIR, "summary.csv"), row.names = FALSE)
  cat("Summary saved to:", file.path(RESULTS_DIR, "summary.csv"), "\n")
}

cat("\nDone.\n")