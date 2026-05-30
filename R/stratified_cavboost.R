###############################################################################
# stratified_cavboost.R
#
# Per-stratum weighted KM gradient for subgroup identification with
# survival outcomes.  Censoring is handled naturally through per-stratum
# risk sets (no IPCW needed — the risk-set structure within each stratum
# already accounts for staggered entry / administrative censoring).
#
# V_strat =  sum_k  W_k   · d_k^(1)  -  sum_k  (n_k − W_k) · d_k^(2)
#
# Gradient for patient j in stratum k (XGBoost minimizes → negated):
#   g_j = −p_j(1−p_j)[ d_k^(1) + d_k^(2) +
#           W_k · ∂d_k^(1)/∂p_j + (n_k−W_k) · ∂d_k^(2)/∂(1-p_j) ]
#
# Reference: Zhang et al. (2020) SubgroupBoost, Stat Med 39(28), 4133–4146.
###############################################################################

library(xgboost)
library(survival)

# =========================================================================
# 1. Per-stratum sorted data
# =========================================================================
make_sorted_stratum <- function(time, status, trt, tau) {
  n <- length(time)
  ord <- order(time)
  list(time = time[ord], status = status[ord], trt = trt[ord],
       ord = ord, iord = order(ord), n = n, tau = tau)
}

# =========================================================================
# 2. Hazard increment and its gradient w.r.t. p
#    w_total is the weight vector (p for subgroup 1, 1−p for subgroup 2).
#    Returns: h = nume / denom,  g = ∂h/∂weight (scalar weight derivative).
#
#    The per-observation gradient ∂h/∂p_j = sign_j · g_j where
#    sign_j = +1 for subgroup 1 (w = p), −1 for subgroup 2 (w = 1−p).
#    The sign is applied by the caller in stratum_gradient_components.
# =========================================================================
hazard_inc_strat <- function(in_e, in_r, w_total, n) {
  denom <- sum(w_total[in_r])
  if (denom <= 0) return(list(h = 0, g = rep(0, n)))
  nume <- sum(w_total[in_e])
  inc <- nume / denom
  # g_j = (I(e_j) − inc · I(r_j)) / denom   (∂h/∂w_j)
  g <- (1 / denom) * (as.numeric(in_e) - inc * as.numeric(in_r))
  list(h = inc, g = g)
}

# =========================================================================
# 3. Compute per-stratum d_k^(1), d_k^(2) and their gradients.
#    All arrays are in the stratum's sorted (by time) order.
# =========================================================================
stratum_gradient_components <- function(p_sorted, sd, tau) {
  n_k <- sd$n
  q_sorted <- 1 - p_sorted
  
  w1 <- p_sorted          # subgroup 1 weights
  w2 <- q_sorted          # subgroup 2 weights
  
  event_times <- unique(sd$time[sd$status == 1])
  event_times <- event_times[event_times <= tau]
  if (length(event_times) == 0) {
    return(list(d1 = 0, g1 = rep(0, n_k), d2 = 0, g2 = rep(0, n_k)))
  }
  
  dt <- c(event_times[1], diff(event_times))
  
  # --- Subgroup 1: p-weighted ---
  H1 <- c(0, 0); gH1 <- list(rep(0, n_k), rep(0, n_k))
  H2 <- c(0, 0); gH2 <- list(rep(0, n_k), rep(0, n_k))
  d1 <- 0; g1 <- rep(0, n_k)
  d2 <- 0; g2 <- rep(0, n_k)
  # Greenwood variance (Aalen form: Var(H(t)) = sum d/Y^2)
  var_H1_trt <- 0; var_H1_ctrl <- 0
  var_H2_trt <- 0; var_H2_ctrl <- 0
  var_d1 <- 0; var_d2 <- 0
  
  for (i in seq_along(event_times)) {
    ti <- event_times[i]; dti <- dt[i]
    in_r1 <- sd$time >= ti & sd$trt == 1
    in_r0 <- sd$time >= ti & sd$trt == 0
    in_e1 <- sd$time == ti & sd$status == 1 & sd$trt == 1
    in_e0 <- sd$time == ti & sd$status == 1 & sd$trt == 0
    
    # Weighted counts for Greenwood
    Y1_trt <- sum(w1[in_r1]); Y1_ctrl <- sum(w1[in_r0])
    d1_trt <- sum(w1[in_e1]); d1_ctrl <- sum(w1[in_e0])
    Y2_trt <- sum(w2[in_r1]); Y2_ctrl <- sum(w2[in_r0])
    d2_trt <- sum(w2[in_e1]); d2_ctrl <- sum(w2[in_e0])
    
    ht <- hazard_inc_strat(in_e1, in_r1, w1, n_k)
    hc <- hazard_inc_strat(in_e0, in_r0, w1, n_k)
    H1 <- H1 + c(ht$h, hc$h)
    gH1[[1]] <- gH1[[1]] + ht$g
    gH1[[2]] <- gH1[[2]] + hc$g
    S1 <- exp(-H1)
    
    ht2 <- hazard_inc_strat(in_e1, in_r1, w2, n_k)
    hc2 <- hazard_inc_strat(in_e0, in_r0, w2, n_k)
    H2 <- H2 + c(ht2$h, hc2$h)
    gH2[[1]] <- gH2[[1]] + ht2$g
    gH2[[2]] <- gH2[[2]] + hc2$g
    S2 <- exp(-H2)
    
    # Accumulate Greenwood (Aalen form: Var(H) += d/Y^2)
    if (Y1_trt > 0) var_H1_trt <- var_H1_trt + d1_trt / (Y1_trt * Y1_trt)
    if (Y1_ctrl > 0) var_H1_ctrl <- var_H1_ctrl + d1_ctrl / (Y1_ctrl * Y1_ctrl)
    if (Y2_trt > 0) var_H2_trt <- var_H2_trt + d2_trt / (Y2_trt * Y2_trt)
    if (Y2_ctrl > 0) var_H2_ctrl <- var_H2_ctrl + d2_ctrl / (Y2_ctrl * Y2_ctrl)
    
    # Var(S(t)) = S(t)^2 * Var(H(t))
    var_d1 <- var_d1 + (S1[1]^2 * var_H1_trt + S1[2]^2 * var_H1_ctrl) * dti^2
    var_d2 <- var_d2 + (S2[1]^2 * var_H2_trt + S2[2]^2 * var_H2_ctrl) * dti^2
    
    d1 <- d1 + (S1[1] - S1[2]) * dti
    g1 <- g1 + (-S1[1] * gH1[[1]] + S1[2] * gH1[[2]]) * dti
    d2 <- d2 + (S2[1] - S2[2]) * dti
    g2 <- g2 + (-S2[1] * gH2[[1]] + S2[2] * gH2[[2]]) * dti
  }
  
  list(d1 = d1, g1 = g1, d2 = d2, g2 = g2, var_d1 = var_d1, var_d2 = var_d2)
}

# =========================================================================
# 4. Custom XGBoost objective (stratified loss)
# =========================================================================
stratified_loss <- function(preds, dtrain) {
  trt <- getinfo(dtrain, "label")
  sorted_strata <- attr(dtrain, "sorted_strata")
  stratum <- attr(dtrain, "stratum")
  K <- length(sorted_strata)
  n <- length(preds)
  
  p_raw <- 1 / (1 + exp(-preds))
  pg <- p_raw * (1 - p_raw)
  
  grad <- numeric(n)
  
  for (k in 1:K) {
    sd <- sorted_strata[[k]]
    if (is.null(sd) || sd$n < 2) next
    
    idx_k <- which(stratum == k)
    n_k <- length(idx_k)
    p_k <- p_raw[idx_k]
    p_sorted <- p_k[sd$ord]
    
    res <- stratum_gradient_components(p_sorted, sd, sd$tau)
    d1 <- res$d1; g1 <- res$g1  # ∂d1/∂p in sorted order
    d2 <- res$d2; g2 <- res$g2  # ∂d2/∂(1-p) in sorted order
    
    W_k <- sum(p_k)
    n_k_minus_W <- n_k - W_k
    
    # NB: + g2 term because ∂d2/∂p = (∂d2/∂(1-p)) · ∂(1-p)/∂p = −g2,
    # and the loss is V_strat = ΣW·d1 − Σ(n−W)·d2, giving ∂V/∂p = d1 + W·g1 + d2 + (n−W)·g2.
    gp_k <- -(d1 + d2 + W_k * g1[sd$iord] + n_k_minus_W * g2[sd$iord])
    grad[idx_k] <- pg[idx_k] * gp_k
  }
  
  list(grad = grad, hess = rep(1.0, n))
}

# =========================================================================
# 4.5 Variance-gated loss: blends original (global) and stratified gradients
#     per stratum using Greenwood variance
# =========================================================================
stratified_loss_vg <- function(preds, dtrain) {
  trt <- getinfo(dtrain, "label")
  sd_orig <- attr(dtrain, "sorted_orig")    # global sorted data (original CAVBoost)
  sorted_strata <- attr(dtrain, "sorted_strata")
  stratum <- attr(dtrain, "stratum")
  K <- length(sorted_strata)
  n <- length(preds)
  
  p_raw <- 1 / (1 + exp(-preds))
  pg <- p_raw * (1 - p_raw)
  
  grad <- numeric(n)
  
  # ---- Global (original) RMST and variance ----
  # Weighted KM on ALL patients (p and 1-p weighted)
  ps <- p_raw[sd_orig$ord]
  et <- unique(sd_orig$time[sd_orig$status == 1])
  et <- et[et <= sd_orig$tau]
  
  d_orig <- c(0, 0); g_orig <- list(rep(0, n), rep(0, n))
  var_H_trt <- 0; var_H_ctrl <- 0; var_H_trt2 <- 0; var_H_ctrl2 <- 0
  var_d_orig <- c(0, 0)
  
  if (length(et) >= 2) {
    dt_et <- c(et[1], diff(et))
    H <- c(0, 0, 0, 0)
    gH <- list(rep(0, n), rep(0, n), rep(0, n), rep(0, n))
    
    for (i in seq_along(et)) {
      ti <- et[i]; dti <- dt_et[i]
      rs <- sd_orig$time >= ti; re <- sd_orig$time == ti & sd_orig$status == 1
      in_r1 <- rs & sd_orig$trt == 1; in_r0 <- rs & sd_orig$trt == 0
      in_e1 <- re & sd_orig$trt == 1; in_e0 <- re & sd_orig$trt == 0
      
      # Weighted counts for Greenwood
      Y_trt <- sum(ps[in_r1]); Y_ctrl <- sum(ps[in_r0])
      d_trt <- sum(ps[in_e1]); d_ctrl <- sum(ps[in_e0])
      Y_trt2 <- sum((1-ps)[in_r1]); Y_ctrl2 <- sum((1-ps)[in_r0])
      d_trt2 <- sum((1-ps)[in_e1]); d_ctrl2 <- sum((1-ps)[in_e0])
      
      ht <- hazard_inc_strat(in_e1, in_r1, ps, n)
      hc <- hazard_inc_strat(in_e0, in_r0, ps, n)
      hh1 <- hazard_inc_strat(in_e1, in_r1, 1 - ps, n)
      hh0 <- hazard_inc_strat(in_e0, in_r0, 1 - ps, n)
      H <- H + c(ht$h, hc$h, hh1$h, hh0$h)
      gH[[1]] <- gH[[1]] + ht$g; gH[[2]] <- gH[[2]] + hc$g
      gH[[3]] <- gH[[3]] + hh1$g; gH[[4]] <- gH[[4]] + hh0$g
      S <- exp(-H)
      
      if (Y_trt > 0) var_H_trt <- var_H_trt + d_trt / (Y_trt * Y_trt)
      if (Y_ctrl > 0) var_H_ctrl <- var_H_ctrl + d_ctrl / (Y_ctrl * Y_ctrl)
      if (Y_trt2 > 0) var_H_trt2 <- var_H_trt2 + d_trt2 / (Y_trt2 * Y_trt2)
      if (Y_ctrl2 > 0) var_H_ctrl2 <- var_H_ctrl2 + d_ctrl2 / (Y_ctrl2 * Y_ctrl2)
      
      var_d_orig[1] <- var_d_orig[1] + (S[1]^2 * var_H_trt + S[2]^2 * var_H_ctrl) * dti^2
      var_d_orig[2] <- var_d_orig[2] + (S[3]^2 * var_H_trt2 + S[4]^2 * var_H_ctrl2) * dti^2
      
      d_orig[1] <- d_orig[1] + (S[1] - S[2]) * dti
      d_orig[2] <- d_orig[2] + (S[3] - S[4]) * dti
      g_orig[[1]] <- g_orig[[1]] + (-S[1] * gH[[1]] + S[2] * gH[[2]]) * dti
      g_orig[[2]] <- g_orig[[2]] + (-S[3] * gH[[3]] + S[4] * gH[[4]]) * dti
    }
  }
  
  sp <- sum(p_raw); sq <- n - sp
  g_orig_vec <- as.numeric(pg * (-(d_orig[1] + sp * g_orig[[1]] + d_orig[2] - sq * g_orig[[2]]))[sd_orig$iord])
  var_orig <- var_d_orig[1] + var_d_orig[2]  # total global variance
  
  # ---- Per-stratum (stratified) gradient and variance ----
  for (k in 1:K) {
    sd <- sorted_strata[[k]]
    if (is.null(sd) || sd$n < 2) next
    
    idx_k <- which(stratum == k)
    n_k <- length(idx_k)
    p_k <- p_raw[idx_k]
    p_sorted <- p_k[sd$ord]
    
    res <- stratum_gradient_components(p_sorted, sd, sd$tau)
    d1 <- res$d1; g1 <- res$g1
    d2 <- res$d2; g2 <- res$g2
    var_d1 <- res$var_d1; var_d2 <- res$var_d2
    
    W_k <- sum(p_k)
    n_k_minus_W <- n_k - W_k
    
    # Stratified gradient
    gp_k <- -(d1 + d2 + W_k * g1[sd$iord] + n_k_minus_W * g2[sd$iord])
    g_strat <- pg[idx_k] * gp_k
    
    # Variance gate: Var(d_k) / (Var(d_orig) + Var(d_k))
    var_k <- var_d1 + var_d2
    alpha_k <- if (var_orig + var_k > 0) var_k / (var_orig + var_k) else 0.5
    alpha_k <- min(max(alpha_k, 0.01), 0.99)  # clamp
    
    # Blend: g = alpha * g_orig + (1 - alpha) * g_strat
    grad[idx_k] <- alpha_k * g_orig_vec[idx_k] + (1 - alpha_k) * g_strat
  }
  
  # Fill in any unassigned indices (non-stratified observations)
  zero_grad <- which(grad == 0 & seq_len(n) > 0)
  unassigned <- setdiff(seq_len(n), unlist(lapply(1:K, function(k) {
    sd <- sorted_strata[[k]]
    if (is.null(sd) || sd$n < 2) return(integer(0))
    which(stratum == k)
  })))
  if (length(unassigned) > 0) {
    grad[unassigned] <- g_orig_vec[unassigned]
  }
  
  list(grad = grad, hess = rep(1.0, n))
}

# =========================================================================
# 5. Prediction
# =========================================================================
pred_stratified <- function(model, dat) {
  features <- setdiff(names(dat), c("trt01p","time","status","A","U","delta","delta_tilde",
                                      "oracle_delta","id","T","C"))
  X <- as.matrix(dat[, features, drop = FALSE])
  1 / (1 + exp(-predict(model, X)))
}

# =========================================================================
# 6. Training
# =========================================================================
# =========================================================================
# 6. Cross-fitted prognostic strata (avoids training-set overlap bias)
# =========================================================================
crossfit_prognostic_strata <- function(dat, features,
                                        time = "time", status = "status",
                                        nfold = 5, K = 4, seed = 42) {
  n <- nrow(dat)
  set.seed(seed)
  folds <- sample(rep(1:nfold, length.out = n))
  lp <- numeric(n)
  
  x <- data.matrix(dat[, features, drop = FALSE])
  y <- survival::Surv(dat[[time]], dat[[status]])
  
  for (fold in 1:nfold) {
    test_idx <- which(folds == fold)
    train_idx <- which(folds != fold)
    cv <- suppressWarnings(glmnet::cv.glmnet(x[train_idx, , drop = FALSE], y[train_idx],
                              family = "cox", alpha = 0.5, nfolds = 5, cox.ties = "breslow"))
    lp[test_idx] <- drop(stats::predict(cv, x[test_idx, , drop = FALSE],
                                          s = "lambda.min"))
  }
  
  qq <- unique(stats::quantile(lp, seq(0, 1, 1 / K), na.rm = TRUE))
  if (length(qq) < 2) return(rep(1, n))
  as.numeric(cut(lp, qq, include.lowest = TRUE, right = TRUE))
}

# =========================================================================
# 7. Training
# =========================================================================
train_stratified_cavboost <- function(dat, time, status, tau, stratum,
                                       eta = 0.05, max_depth = 4, nr = 50,
                                       nthread = 1L) {
  if ("A" %in% names(dat) && !("trt01p" %in% names(dat))) dat$trt01p <- dat$A
  if ("U" %in% names(dat) && !("time" %in% names(dat))) time <- dat$U
  if ("delta_tilde" %in% names(dat) && !("status" %in% names(dat))) status <- dat$delta_tilde
  
  # Pre-sort data per stratum (no IPCW — per-stratum risk sets handle censoring)
  K <- max(stratum)
  sorted_strata <- vector("list", K)
  for (k in 1:K) {
    idx <- which(stratum == k)
    if (length(idx) < 2) { sorted_strata[[k]] <- NULL; next }
    sorted_strata[[k]] <- make_sorted_stratum(
      time[idx], status[idx], dat$trt01p[idx], tau)
  }
  
  features <- setdiff(names(dat), c("trt01p","time","status","A","U","delta","delta_tilde",
                                      "oracle_delta","id","T","C"))
  Xmat <- as.matrix(dat[, features, drop = FALSE])
  
  dtrain <- xgb.DMatrix(Xmat, label = dat$trt01p)
  attr(dtrain, "sorted_strata") <- sorted_strata
  attr(dtrain, "stratum") <- stratum
  
  params <- list(eta = eta, max_depth = max_depth, lambda = 1,
                 min_child_weight = 0, subsample = 1, colsample_bytree = 1,
                 nthread = nthread)
  
  xgb.train(params = params, data = dtrain, nrounds = nr,
             objective = stratified_loss, verbose = 0)
}

# =========================================================================
# 8. Training with variance-gated gradient (blends orig + stratified)
# =========================================================================
train_stratified_cavboost_vg <- function(dat, time, status, tau, stratum,
                                          eta = 0.1, max_depth = 2, nr = 50) {
  if ("A" %in% names(dat) && !("trt01p" %in% names(dat))) dat$trt01p <- dat$A
  if ("U" %in% names(dat) && !("time" %in% names(dat))) time <- dat$U
  if ("delta_tilde" %in% names(dat) && !("status" %in% names(dat))) status <- dat$delta_tilde
  
  # Global sorted data (original CAVBoost style, no IPCW)
  sd_orig <- make_sorted_data(time, status, dat$trt01p, tau)
  
  # Per-stratum sorted data
  K <- max(stratum)
  sorted_strata <- vector("list", K)
  for (k in 1:K) {
    idx <- which(stratum == k)
    if (length(idx) < 2) { sorted_strata[[k]] <- NULL; next }
    sorted_strata[[k]] <- make_sorted_stratum(
      time[idx], status[idx], dat$trt01p[idx], tau)
  }
  
  features <- setdiff(names(dat), c("trt01p","time","status","A","U","delta","delta_tilde",
                                      "oracle_delta","id","T","C"))
  Xmat <- as.matrix(dat[, features, drop = FALSE])
  
  dtrain <- xgb.DMatrix(Xmat, label = dat$trt01p)
  attr(dtrain, "sorted_orig") <- sd_orig
  attr(dtrain, "sorted_strata") <- sorted_strata
  attr(dtrain, "stratum") <- stratum
  
  params <- list(eta = eta, max_depth = max_depth, lambda = 1,
                 min_child_weight = 0, subsample = 1, colsample_bytree = 1)
  
  xgb.train(params = params, data = dtrain, nrounds = nr,
             objective = stratified_loss_vg, verbose = 0)
}

# =========================================================================
# 9. Holdout model selection: train both original and stratified,
#    select the one with better held-out RMST gain, refit on full data.
# =========================================================================

# Helper: compute realized RMST gain from a binary subgroup rule
rmst_gain <- function(pred, time, status, tau) {
  sub_g <- which(pred > 0.5)
  comp <- which(pred <= 0.5)
  if (length(sub_g) < 3 || length(comp) < 3) return(0)
  
  sf <- survfit(Surv(time[sub_g], status[sub_g]) ~ 1)
  ss <- summary(sf, times = seq(0, tau, length.out = 200), extend = TRUE)
  rmst_sub <- mean(ss$surv, na.rm = TRUE) * tau
  
  sf <- survfit(Surv(time[comp], status[comp]) ~ 1)
  ss <- summary(sf, times = seq(0, tau, length.out = 200), extend = TRUE)
  rmst_comp <- mean(ss$surv, na.rm = TRUE) * tau
  
  rmst_sub - rmst_comp
}

select_model_by_holdout <- function(dat, time, status, tau, stratum, features,
                                     params_orig = list(eta = 0.05, max_depth = 3, nr = 50),
                                     params_strat = list(eta = 0.1, max_depth = 2, nr = 50)) {
  n <- nrow(dat)
  set.seed(42)
  train_idx <- sample(1:n, round(0.7 * n))
  hold_idx <- setdiff(1:n, train_idx)
  
  dt_tr <- dat[train_idx, ]
  dt_hold <- dat[hold_idx, ]
  
  # Train original on 70%
  fo <- tryCatch(
    train_rmst_cavboost(dt_tr, dt_tr[[time]], dt_tr[[status]], tau,
                         eta = params_orig$eta, max_depth = params_orig$max_depth,
                         nr = params_orig$nr, covars = NULL),
    error = function(e) NULL)
  
  # Train stratified on 70%
  st_hold <- stratum[train_idx]
  fs <- tryCatch(
    train_stratified_cavboost(dt_tr, dt_tr[[time]], dt_tr[[status]], tau,
                               stratum = st_hold,
                               eta = params_strat$eta,
                               max_depth = params_strat$max_depth,
                               nr = params_strat$nr),
    error = function(e) NULL)
  
  # Evaluate on holdout
  gain_orig <- gain_strat <- -Inf
  
  if (!is.null(fo)) {
    po <- pred_subgroup(fo, dt_hold)
    gain_orig <- rmst_gain(po, dt_hold[[time]], dt_hold[[status]], tau)
  }
  if (!is.null(fs)) {
    ps <- pred_stratified(fs, dt_hold)
    gain_strat <- rmst_gain(ps, dt_hold[[time]], dt_hold[[status]], tau)
  }
  
  # Pick winner
  winner <- if (gain_orig >= gain_strat) "original" else "stratified"
  
  # Refit on full data
  if (winner == "original") {
    fit <- train_rmst_cavboost(dat, dat[[time]], dat[[status]], tau,
                                eta = params_orig$eta, max_depth = params_orig$max_depth,
                                nr = params_orig$nr, covars = NULL)
  } else {
    fit <- train_stratified_cavboost(dat, dat[[time]], dat[[status]], tau,
                                      stratum = stratum,
                                      eta = params_strat$eta,
                                      max_depth = params_strat$max_depth,
                                      nr = params_strat$nr)
  }
  
  list(fit = fit, winner = winner,
       gain_orig = gain_orig, gain_strat = gain_strat)
}
