# Main comparison: Orig vs Strat_CF vs VT, 50 reps, 4 scenarios
# Fresh R process per rep, parallel per scenario
# Usage: Rscript run_main_comparison_50.R

library(pROC); library(mvtnorm)
source("R/stratified_cavboost.R")
source("~/.openclaw/workspace/CAVBoost/rmst_cavboost_clean.R")

Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT = "600")

tau <- 30; n_tr <- 500; n_te <- 2000; rho <- 1/3; nr <- 50
Smat <- matrix(rho, 52, 52); diag(Smat) <- 1


run_one <- function(sc, rep) {
  set.seed(42 + rep * 1000 + sc * 100)
  X <- rmvnorm(n_tr + n_te, sigma = Smat)
  colnames(X) <- c(paste0("z", 1:50), "S1", "S2")
  colnames(X)[1:4] <- c("Z1","Z2","Z3","Z4")
  A <- rbinom(n_tr + n_te, 1, 0.5); b0 <- sqrt(6); s0 <- 0.4
  Zb <- X[, 1:4, drop = FALSE] %*% c(0.4, 0.4, 0.4, 0.4)
  te <- switch(sc,
    `1` = X[, "S1"],
    `2` = X[, "S1"] - X[, "S2"],
    `3` = { s <- X[, "S1"]; 2 * ifelse(abs(s) < 0.67, exp(-s^2) - 0.4, exp(-s^2) - 0.8) },
    `4` = { s1 <- X[, "S1"]; 2 * ((-1.07 <= s1 & s1 < 1.07) & (-1.07 <= X[, "S2"] & X[, "S2"] < 1.07)) - 1 },
    `5` = { s <- X[, "S1"]; 2 * ifelse(s >= 0.67 | (-0.67 <= s & s < 0), 1, 0) - 1 },
    `6` = { s1 <- X[, "S1"]; s2 <- X[, "S2"]; 2 * ifelse((s1 >= 0 & s2 >= -0.67) | (s1 < 0 & s2 < -0.67), 1, 0) - 1 }
  )
  if (sc %in% 4:6) Zb <- -Zb^2
  T <- exp(b0 + A * te + Zb + s0 * rnorm(n_tr + n_te))
  C <- pmin(30, rexp(n_tr + n_te, rate = -log(0.9) / 12))
  U <- pmin(T, C, tau); st <- as.numeric(T <= C)
  oracle <- as.numeric(te > 0)
  if (length(unique(oracle)) < 2) return(NULL)
  
  tr <- data.frame(X[1:n_tr, ], trt01p = A[1:n_tr], time = U[1:n_tr], status = st[1:n_tr])
  te_df <- data.frame(X[-(1:n_tr), ], trt01p = A[-(1:n_tr)], time = U[-(1:n_tr)], status = st[-(1:n_tr)])
  l <- oracle[-(1:n_tr)]
  
  auc_ <- function(p, l) tryCatch(as.numeric(pROC::auc(pROC::roc(l, p, quiet = TRUE, direction = "<"))), error = function(e) NA)
  
  # Original
  fo <- tryCatch(train_rmst_cavboost(tr, tr$time, tr$status, tau, eta = 0.05, max_depth = 3, nr = nr), error = function(e) NULL)
  po <- if (!is.null(fo)) pred_subgroup(fo, te_df) else rep(0.5, nrow(te_df))
  rm(fo)
  
  # Cross-fitted stratified (6 prognostic vars)
  feat_all <- colnames(X)  # all 52 covariates
  st_ <- tryCatch(crossfit_prognostic_strata(tr, feat_all, nfold = 5, K = 4, seed = rep * 100 + sc), error = function(e) NULL)
  fs <- if (!is.null(st_)) tryCatch(train_stratified_cavboost(tr, tr$time, tr$status, tau, stratum = st_, eta = 0.1, max_depth = 2, nr = nr), error = function(e) NULL) else NULL
  ps <- if (!is.null(fs)) pred_stratified(fs, te_df) else rep(0.5, nrow(te_df))
  rm(fs)
  
  # VT: Cox per arm
  fm <- paste("Surv(time,status)~", paste(colnames(X), collapse = "+"))
  ctrl <- tr[tr$trt01p == 0, ]; trt_d <- tr[tr$trt01p == 1, ]
  fc <- tryCatch(coxph(as.formula(fm), data = ctrl, x = TRUE), error = function(e) NULL)
  ft <- tryCatch(coxph(as.formula(fm), data = trt_d, x = TRUE), error = function(e) NULL)
  vt <- if (!is.null(fc) && !is.null(ft)) predict(ft, te_df, type = "lp") - predict(fc, te_df, type = "lp") else rep(0, nrow(te_df))
  
  data.frame(
    sc = sc, rep = rep,
    orig_auc = auc_(po, l), strat_auc = auc_(ps, l), vt_auc = auc_(vt, l)
  )
}

# Run
sc_names <- list("1" = "S1_Linear", "2" = "S2_Diff", "3" = "S3_U", "4" = "S4_Enclave", "5" = "S5_S", "6" = "S6_Cross")
all <- data.frame()
n_sim <- 50

for (sc in c(1, 2, 3, 4, 5, 6)) {
  res <- data.frame()
  for (rep in 1:n_sim) {
    r <- tryCatch(run_one(sc, rep), error = function(e) NULL)
    if (!is.null(r)) res <- rbind(res, r)
    cat(sprintf("\r%s rep %d/%d", sc_names[[as.character(sc)]], rep, n_sim))
  }
  cat("\n")
  all <- rbind(all, res)
  saveRDS(all, "main_comparison_50rep.rds")
}

cat("\n============================================================\n")
cat("  MAIN COMPARISON (50 reps, cross-fitted)\n")
cat("============================================================\n\n")
for (sc_name in c("S1_Linear", "S2_Diff", "S3_U", "S4_Enclave", "S5_S", "S6_Cross")) {
  s <- all[all$sc == sc_name, ]
  cat(sprintf("--- %s (n=%d) ---\n", sc_name, nrow(s)))
  cat(sprintf("%-10s %8s\n", "Method", "AUC"))
  for (mn in list(c("orig", "Original"), c("strat", "Strat_CF"), c("vt", "VT"))) {
    cat(sprintf("%-10s %8.4f\n", mn[[2]], mean(s[[paste0(mn[[1]], "_auc")]], na.rm = TRUE)))
  }
  cat("\n")
}
cat("Saved to main_comparison_50rep.rds\n")
