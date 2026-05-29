# Variance-gated gradient: blends original + stratified per stratum
# Gate: alpha_k = Var(d_k) / (Var(d_orig) + Var(d_k))
# g_j = alpha_k * g_orig_j + (1 - alpha_k) * g_strat_j
# Usage: Rscript run_variance_gate.R  (from repo root)

set.seed(42)
library(survival); library(xgboost); library(mvtnorm); library(pROC)
source("R/stratified_cavboost.R")
source("~/.openclaw/workspace/CAVBoost/rmst_cavboost_clean.R")

tau <- 30; n_sim <- 20; n_tr <- 500; n_te <- 2000; rho <- 1/3; nr <- 30
Smat <- matrix(rho, 52, 52); diag(Smat) <- 1

run_vg <- function(scenario, sc) {
  r_o <- r_s <- r_vg <- list()
  for (rep in 1:n_sim) {
    set.seed(42 + rep * 1000 + sc * 100)
    X <- rmvnorm(n_tr + n_te, sigma = Smat)
    colnames(X) <- c(paste0("z", 1:50), "S1", "S2")
    colnames(X)[1:4] <- c("Z1","Z2","Z3","Z4")
    A <- rbinom(n_tr + n_te, 1, 0.5); b0 <- sqrt(6); s0 <- 0.4
    Zb <- X[, 1:4, drop = FALSE] %*% c(0.4, 0.4, 0.4, 0.4)
    te <- if (sc == 4) {
      s1 <- X[, "S1"]
      2 * ((-1.07 <= s1 & s1 < 1.07) & (-1.07 <= X[, "S2"] & X[, "S2"] < 1.07)) - 1
    } else {
      s <- X[, "S1"]
      2 * ifelse(s >= 0.67 | (-0.67 <= s & s < 0), 1, 0) - 1
    }
    Zb <- -Zb^2
    T <- exp(b0 + A * te + Zb + s0 * rnorm(n_tr + n_te))
    C <- pmin(30, rexp(n_tr + n_te, rate = -log(0.9) / 12))
    U <- pmin(T, C, tau); st <- as.numeric(T <= C); oracle <- as.numeric(te > 0)
    if (length(unique(oracle)) < 2) next

    tr <- data.frame(X[1:n_tr, ], trt01p = A[1:n_tr], time = U[1:n_tr], status = st[1:n_tr])
    te_df <- data.frame(X[-(1:n_tr), ], trt01p = A[-(1:n_tr)], time = U[-(1:n_tr)], status = st[-(1:n_tr)])
    l <- oracle[-(1:n_tr)]; f <- colnames(X)
    feat_fit <- c("Z1","Z2","Z3","Z4","S1","S2")

    # Original CAVBoost
    fo <- tryCatch(train_rmst_cavboost(tr, tr$time, tr$status, tau, eta=0.05, max_depth=3, nr=nr, covars=NULL),
                    error = function(e) NULL)
    if (is.null(fo)) { cat(sprintf("orig fail %s %d\n", scenario, rep)); next }
    po <- pred_subgroup(fo, te_df)

    # Cross-fitted stratified (no gate)
    st_ <- tryCatch(crossfit_prognostic_strata(tr, feat_fit, nfold=5, K=4, seed=rep*100+sc),
                     error = function(e) NULL)
    fs <- if (!is.null(st_)) tryCatch(train_stratified_cavboost(tr, tr$time, tr$status, tau, stratum = st_,
                                                eta = 0.1, max_depth = 2, nr = nr),
                      error = function(e) NULL) else NULL
    ps_te <- if (!is.null(fs)) pred_stratified(fs, te_df) else rep(0.5, nrow(te_df))

    # Variance-gated training
    fvg <- if (!is.null(st_)) tryCatch(train_stratified_cavboost_vg(tr, tr$time, tr$status, tau, stratum = st_,
                                                eta = 0.1, max_depth = 2, nr = nr),
                      error = function(e) NULL) else NULL
    pvg_te <- if (!is.null(fvg)) pred_stratified(fvg, te_df) else rep(0.5, nrow(te_df))

    cmf <- function(p, ll) {
      r <- suppressMessages(roc(ll, p, quiet = TRUE))
      c(AUC = as.numeric(auc(r)), Acc = mean((p > .5) == ll),
        FPR = if (sum(!ll) > 0) sum(p > .5 & !ll) / sum(!ll) else NA,
        FNR = if (sum(ll) > 0) sum(p <= .5 & ll) / sum(ll) else NA)
    }
    r_o[[rep]] <- cmf(po, l)
    r_s[[rep]] <- cmf(ps_te, l)
    r_vg[[rep]] <- cmf(pvg_te, l)
  }
  to_mat <- function(lst) {
    ok <- !sapply(lst, is.null); if (sum(ok) == 0) return(NULL); do.call(rbind, lst[ok])
  }
  ro <- to_mat(r_o); rs <- to_mat(r_s); rvg <- to_mat(r_vg)
  rn <- function(d) {
    if (is.null(d)) return(list(NA_real_, NA_real_, NA_real_, NA_real_, 0L))
    list(mean(d[, 1], na.rm = TRUE), mean(d[, 2], na.rm = TRUE),
         mean(d[, 3], na.rm = TRUE), mean(d[, 4], na.rm = TRUE), nrow(d))
  }
  cat(sprintf("\n=== %s (var-gate, %d reps) ===\n", scenario, n_sim))
  cat(sprintf("%-15s %8s %8s %8s %8s\n", "Method", "AUC", "Acc", "FPR", "FNR"))
  for (lst in list(list("Orig", ro), list("Strat_CF", rs), list("VarGate", rvg))) {
    nm <- lst[[1]]; d <- lst[[2]]; rr <- rn(d)
    cat(sprintf("%-15s %8.4f %8.4f %8.4f %8.4f (N=%d)\n",
                nm, rr[[1]], rr[[2]], rr[[3]], rr[[4]], rr[[5]]))
  }
}
run_vg("S4_Enclave", 4)
run_vg("S5_SShaped", 5)
