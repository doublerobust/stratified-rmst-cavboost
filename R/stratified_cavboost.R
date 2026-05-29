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
hazard_inc <- function(in_e, in_r, w_total, n) {
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
  d1 <- 0; g1 <- rep(0, n_k)
  for (i in seq_along(event_times)) {
    ti <- event_times[i]; dti <- dt[i]
    in_r1 <- sd$time >= ti & sd$trt == 1
    in_r0 <- sd$time >= ti & sd$trt == 0
    in_e1 <- sd$time == ti & sd$status == 1 & sd$trt == 1
    in_e0 <- sd$time == ti & sd$status == 1 & sd$trt == 0
    ht <- hazard_inc(in_e1, in_r1, w1, n_k)
    hc <- hazard_inc(in_e0, in_r0, w1, n_k)
    H1 <- H1 + c(ht$h, hc$h)
    gH1[[1]] <- gH1[[1]] + ht$g
    gH1[[2]] <- gH1[[2]] + hc$g
    S <- exp(-H1)
    d1 <- d1 + (S[1] - S[2]) * dti
    g1 <- g1 + (-S[1] * gH1[[1]] + S[2] * gH1[[2]]) * dti
  }
  
  # --- Subgroup 2: (1-p)-weighted ---
  H2 <- c(0, 0); gH2 <- list(rep(0, n_k), rep(0, n_k))
  d2 <- 0; g2 <- rep(0, n_k)
  for (i in seq_along(event_times)) {
    ti <- event_times[i]; dti <- dt[i]
    in_r1 <- sd$time >= ti & sd$trt == 1
    in_r0 <- sd$time >= ti & sd$trt == 0
    in_e1 <- sd$time == ti & sd$status == 1 & sd$trt == 1
    in_e0 <- sd$time == ti & sd$status == 1 & sd$trt == 0
    ht <- hazard_inc(in_e1, in_r1, w2, n_k)
    hc <- hazard_inc(in_e0, in_r0, w2, n_k)
    H2 <- H2 + c(ht$h, hc$h)
    gH2[[1]] <- gH2[[1]] + ht$g
    gH2[[2]] <- gH2[[2]] + hc$g
    S <- exp(-H2)
    d2 <- d2 + (S[1] - S[2]) * dti
    g2 <- g2 + (-S[1] * gH2[[1]] + S[2] * gH2[[2]]) * dti
  }
  
  list(d1 = d1, g1 = g1, d2 = d2, g2 = g2)
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
train_stratified_cavboost <- function(dat, time, status, tau, stratum,
                                       eta = 0.05, max_depth = 4, nr = 50) {
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
                 min_child_weight = 0, subsample = 1, colsample_bytree = 1)
  
  xgb.train(params = params, data = dtrain, nrounds = nr,
             objective = stratified_loss, verbose = 0)
}
