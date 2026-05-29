###############################################################################
# Alternating gradient test: original CAVBoost ↔ stratified (per round)
# Uses 6 prognostic variables (Z1-Z4, S1, S2) for the Cox PH score
###############################################################################
set.seed(42)
library(survival); library(xgboost); library(mvtnorm)
library(pROC)
source("R/stratified_cavboost.R")
source("../CAVBoost/rmst_cavboost_clean.R")
tau <- 30; n_sim <- 20; n_tr <- 500; n_te <- 2000; rho <- 1/3; nr <- 30
Smat <- matrix(rho, 52, 52); diag(Smat) <- 1; prog_vars <- c("Z1","Z2","Z3","Z4","S1","S2")

run_alt <- function(scenario, sc) {
  r_o <- r_s <- r_a <- list()
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

    # Original CAVBoost
    fo <- tryCatch(train_rmst_cavboost(tr, tr$time, tr$status, tau, eta = 0.05, max_depth = 3, nr = nr, covars = NULL),
                    error = function(e) NULL)
    if (is.null(fo)) { cat(sprintf("orig fail %s %d\n", scenario, rep)); next }
    po <- pred_subgroup(fo, te_df)

    # Stratified (cross-fitted Cox on 6 vars only, not all 52)
    feat_fit <- c("Z1","Z2","Z3","Z4","S1","S2")
    st_ <- tryCatch(crossfit_prognostic_strata(tr, feat_fit, nfold=5, K=4, seed=rep*100+sc),
                     error = function(e) NULL)
    fs <- if (!is.null(st_)) tryCatch(train_stratified_cavboost(tr, tr$time, tr$status, tau, stratum = st_,
                                                eta = 0.1, max_depth = 2, nr = nr),
                      error = function(e) NULL) else NULL
    ps_te <- if (!is.null(fs)) pred_stratified(fs, te_df) else rep(0.5, nrow(te_df))

    # Alternating: build custom eta vector round by round (no IPCW)
    sd_orig <- make_sorted_data(tr$time, tr$status, rep(1, n_tr), tr$trt01p, tau)
    K <- 4
    sorted_strata <- vector("list", K)
    for (k in 1:K) {
      idx <- which(st_ == k)
      if (length(idx) >= 2) sorted_strata[[k]] <- make_sorted_stratum(tr$time[idx], tr$status[idx], tr$trt01p[idx], tau)
    }
    eta <- rep(0, n_tr); tree_list <- list()
    for (rd in 1:nr) {
      p <- 1/(1 + exp(-eta)); pg <- p * (1 - p)
      if (rd %% 2 == 1) {
        # Original gradient (global p-weighted KM, no IPCW)
        ps <- p[sd_orig$ord]
        et <- unique(sd_orig$time[sd_orig$status == 1]); et <- et[et <= tau]
        if (length(et) < 2) { g <- rep(0, n_tr) } else {
          dt_et <- c(et[1], diff(et)); H <- c(0, 0, 0, 0)
          gH <- list(rep(0, n_tr), rep(0, n_tr), rep(0, n_tr), rep(0, n_tr))
          d <- c(0, 0); g <- list(rep(0, n_tr), rep(0, n_tr))
          for (i in seq_along(et)) {
            ti <- et[i]; dti <- dt_et[i]; rs <- sd_orig$time >= ti; re <- sd_orig$time == ti & sd_orig$status == 1
            in_r1 <- rs & sd_orig$trt == 1; in_r0 <- rs & sd_orig$trt == 0
            in_e1 <- re & sd_orig$trt == 1; in_e0 <- re & sd_orig$trt == 0
            ht <- hazard_inc(in_e1, in_r1, ps, n_tr)
            hc <- hazard_inc(in_e0, in_r0, ps, n_tr)
            hh1 <- hazard_inc(in_e1, in_r1, 1 - ps, n_tr)
            hh0 <- hazard_inc(in_e0, in_r0, 1 - ps, n_tr)
            H <- H + c(ht$h, hc$h, hh1$h, hh0$h)
            gH[[1]] <- gH[[1]] + ht$g; gH[[2]] <- gH[[2]] + hc$g
            gH[[3]] <- gH[[3]] + hh1$g; gH[[4]] <- gH[[4]] + hh0$g
            S <- exp(-H); d[1] <- d[1] + (S[1] - S[2]) * dti; d[2] <- d[2] + (S[3] - S[4]) * dti
            g[[1]] <- g[[1]] + (-S[1] * gH[[1]] + S[2] * gH[[2]]) * dti
            g[[2]] <- g[[2]] + (-S[3] * gH[[3]] + S[4] * gH[[4]]) * dti
          }
          sp <- sum(p); sq <- sum(1 - p)
          g <- as.numeric(pg * (-(d[1] + sp * g[[1]] + d[2] - sq * g[[2]]))[sd_orig$iord])
        }
      } else {
        # Stratified gradient (no IPCW)
        g <- rep(0, n_tr)
        for (k in 1:K) {
          sd <- sorted_strata[[k]]
          if (is.null(sd) || sd$n < 2) next
          idx_k <- which(st_ == k)
          pk <- p[idx_k]
          rr <- stratum_gradient_components(pk[sd$ord], sd, tau)
          Wk <- sum(pk)
          # Correct sign: + (n_k - Wk) * g2 (see stratified_loss in stratified_cavboost.R)
          g[idx_k] <- pg[idx_k] * (-(rr$d1 + rr$d2 + Wk * rr$g1[sd$iord] + (length(pk) - Wk) * rr$g2[sd$iord]))
        }
      }
      dt <- xgb.DMatrix(as.matrix(tr[, f]), label = g)
      bt <- xgb.train(params = list(eta = 0.05, max_depth = 2), data = dt, nrounds = 1,
                       verbose = 0, objective = "reg:squarederror")
      eta <- eta + predict(bt, as.matrix(tr[, f]))
      tree_list[[rd]] <- bt
    }
    # Predict on test: sum tree outputs
    pa_te <- 1/(1 + exp(-rowSums(sapply(tree_list, function(b) predict(b, as.matrix(te_df[, f]))))))

    cmf <- function(p, ll) {
      r <- roc(ll, p, direction = "auto")
      c(AUC = as.numeric(r$auc), Acc = mean((p > .5) == ll),
        FPR = if (sum(!ll) > 0) sum(p > .5 & !ll) / sum(!ll) else NA,
        FNR = if (sum(ll) > 0) sum(p <= .5 & ll) / sum(ll) else NA)
    }
    r_o[[rep]] <- cmf(po, l); r_s[[rep]] <- cmf(ps_te, l); r_a[[rep]] <- cmf(pa_te, l)
  }
  # Collapse with NULL-safe rbind
  to_mat <- function(lst) {
    ok <- !sapply(lst, is.null)
    if (sum(ok) == 0) return(NULL)
    do.call(rbind, lst[ok])
  }
  ro <- to_mat(r_o); rs <- to_mat(r_s); ra <- to_mat(r_a)
  rn <- function(d) {
    if (is.null(d)) return(list(NA_real_, NA_real_, NA_real_, NA_real_, 0L))
    list(mean(d[, 1], na.rm = TRUE), mean(d[, 2], na.rm = TRUE),
         mean(d[, 3], na.rm = TRUE), mean(d[, 4], na.rm = TRUE), nrow(d))
  }
  cat(sprintf("\n=== %s (alternating, %d reps) ===\n", scenario, n_sim))
  cat(sprintf("%-15s %8s %8s %8s %8s\n", "Method", "AUC", "Acc", "FPR", "FNR"))
  for (lst in list(c("Orig", "ro"), c("Strat", "rs"), c("Alternating", "ra"))) {
    nm <- lst[1]; d <- get(lst[2]); rr <- rn(d)
    cat(sprintf("%-15s %8.4f %8.4f %8.4f %8.4f (N=%d)\n",
                nm, rr[[1]], rr[[2]], rr[[3]], rr[[4]], rr[[5]]))
  }
}
run_alt("S4_Enclave", 4)
run_alt("S5_SShaped", 5)
