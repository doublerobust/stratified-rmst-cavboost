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
  library(ranger)
})

# ---- Config (defaults; override via N_CONFIGS <<- before source()) ----
if (!exists("N_CONFIGS")) N_CONFIGS <- 1000
if (!exists("N_REPS")) N_REPS <- 5
if (!exists("N_WORKERS")) N_WORKERS <- 11
if (!exists("PARALLEL")) PARALLEL <- TRUE
N_TOTAL <- N_CONFIGS * N_REPS
TAU <- 30
SEED_BASE <- 20260601
NR <- 30  # boosting iterations
NTHREAD <- 2  # XGBoost threads per worker (N_WORKERS * NTHREAD <= total cores)

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
      te_scale = cfg$te_scale,
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

    # ---- Fit VT (ranger per arm) ----
    ctrl <- train[train$trt01p == 0, ]
    trt_d <- train[train$trt01p == 1, ]
    vt_auc <- NA_real_
    vt_preds <- NULL
    if (nrow(ctrl) >= 20 && nrow(trt_d) >= 20) {
      zcols_vt <- setdiff(names(train), c("trt01p", "time", "status"))
      rf_c <- tryCatch(
        ranger(Surv(time, status) ~ ., data = ctrl[, c("time", "status", zcols_vt)],
               num.trees = 200, min.node.size = 10, seed = seed + 300),
        error = function(e) NULL
      )
      rf_t <- tryCatch(
        ranger(Surv(time, status) ~ ., data = trt_d[, c("time", "status", zcols_vt)],
               num.trees = 200, min.node.size = 10, seed = seed + 301),
        error = function(e) NULL
      )
      if (!is.null(rf_c) && !is.null(rf_t)) {
        pc <- predict(rf_c, test[, zcols_vt])
        pt <- predict(rf_t, test[, zcols_vt])
        # RMST via mean survival on fine grid
        tg <- seq(0, tau, length.out = 200)
        surv_grid <- function(surv_mat, times, grid) {
          apply(surv_mat, 1, function(s) {
            stats::approx(times, s, grid, rule = 2, yleft = 1)$y
          })
        }
        sc <- surv_grid(pc$survival, pc$unique.death.times, tg)
        st <- surv_grid(pt$survival, pt$unique.death.times, tg)
        vt_preds <- colMeans(st) * tau - colMeans(sc) * tau
        vt_auc <- auc_(vt_preds, oracle)
      }
    }

    # ---- Build results ----
    aucs_all <- c(aucs, vt_auc)
    names(aucs_all) <- c(paste0("auc_K", k_values), "auc_VT")

    result <- list(
      config_idx = config_idx,
      config = cfg,
      rep = rep,
      seed = seed,
      aucs = aucs_all,
      oracle_optimal_method = which.max(aucs_all),
      vt_auc = vt_auc,
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

# Determine config range (used for multi-process: each process handles a slice)
config_indices <- seq_len(N_CONFIGS)
if (exists("CONFIG_START") && !is.na(CONFIG_START) &&
    exists("CONFIG_END") && !is.na(CONFIG_END)) {
  config_indices <- config_indices[config_indices >= CONFIG_START & config_indices <= CONFIG_END]
  cat(sprintf("Processing configs %d .. %d (%d configs)\n\n",
              CONFIG_START, CONFIG_END, length(config_indices)))
} else {
  cat(sprintf("Processing all %d configs sequentially\n\n", N_CONFIGS))
}
utils::flush.console()

lapply(config_indices, process_rep,
       all_configs = configs,
       moe_dir = MOE_DIR,
       repo_dir = REPO_DIR,
       raw_dir = RAW_DIR,
       results_dir = RESULTS_DIR,
       tau = TAU,
       nr = NR,
       nthread = 1L)

cat("\n=== Simulation Complete ===\n")

# ---- Summarize ----
result_files <- list.files(RESULTS_DIR, "rep_.*\\.rds$", full.names = TRUE)
cat(sprintf("Files written: %d / %d\n", length(result_files), N_TOTAL))

if (length(result_files) > 0) {
  summary_list <- lapply(result_files, function(f) {
    r <- readRDS(f)
    aucs <- r$aucs
    data.frame(
      family = r$config$family,
      n_train = r$config$n_train,
      auc_K1 = if (length(aucs) >= 1) aucs[1] else NA,
      auc_K2 = if (length(aucs) >= 2) aucs[2] else NA,
      auc_K3 = if (length(aucs) >= 3) aucs[3] else NA,
      auc_K4 = if (length(aucs) >= 4) aucs[4] else NA,
      auc_K5 = if (length(aucs) >= 5) aucs[5] else NA,
      auc_VT = if (length(aucs) >= 6) aucs[6] else NA,
      oracle_optimal_method = if (!is.null(r$oracle_optimal_method)) r$oracle_optimal_method
                          else if (!is.null(r$oracle_optimal_K)) r$oracle_optimal_K else NA,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_list)

  cat("\nOptimal method distribution (overall):\n")
  print(table(summary_df$oracle_optimal_method))

  cat("\nVT wins when? (oracle_optimal_method == 6)\n")
  vt_wins <- summary_df$oracle_optimal_method == 6
  cat(sprintf("VT wins: %d / %d (%.1f%%)\n", sum(vt_wins, na.rm = TRUE), nrow(summary_df),
      100 * mean(vt_wins, na.rm = TRUE)))

  cat("\nOptimal method by family:\n")
  method_labels <- c("K1","K2","K3","K4","K5","VT")
  tbl <- table(summary_df$family, summary_df$oracle_optimal_method)
  colnames(tbl) <- method_labels[as.numeric(colnames(tbl))]
  print(round(prop.table(tbl, 1) * 100, 1))

  cat("\nOptimal method by n_train:\n")
  tbl_n <- table(summary_df$n_train, summary_df$oracle_optimal_method)
  colnames(tbl_n) <- method_labels[as.numeric(colnames(tbl_n))]
  print(round(prop.table(tbl_n, 1) * 100, 1))

  # VT AUC stats
  vt_aucs <- summary_df$auc_VT
  cat(sprintf("\nVT AUC: mean = %.4f, median = %.4f, %% missing = %.1f%%\n",
      mean(vt_aucs, na.rm = TRUE), median(vt_aucs, na.rm = TRUE),
      100 * mean(is.na(vt_aucs))))

  # Flag near-perfect AUCs for diagnostic
  auc_cols <- grep("^auc_", names(summary_df), value = TRUE)
  perfect <- apply(summary_df[, auc_cols] > 0.999, 1, any, na.rm = TRUE)
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