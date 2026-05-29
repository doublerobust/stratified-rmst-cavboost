###############################################################################
# stratified_cavboost.R (CORRECTED)
# Per-stratum weighted KM gradient, matching original CAVBoost structure.
#
# V_strat = Σ_k W_k · d_k^(1) − Σ_k (n_k − W_k) · d_k^(2)
#
# Gradient for patient j in stratum k (XGBoost minimizes → negated):
#   g_j = −p_j(1−p_j)[ d_k^(1) + d_k^(2) +
#           W_k · ∂d_k^(1)/∂p_j − (n_k−W_k) · ∂d_k^(2)/∂p_j ]
###############################################################################

library(xgboost)
library(survival)

# =========================================================================
# 1. IPCW weights (same as rmst_cavboost_clean.R)
# =========================================================================
compute_ipcw_weights <- function(time, status, tau, covars = NULL, trt = NULL, eps = 0.05) {
  Y <- pmin(time, tau)
  delta_tilde <- status | (time >= tau)
  cens_ind <- 1 - status
  if (is.null(covars)) {
    fit_cens <- survfit(Surv(time, cens_ind) ~ 1)
    G_hat <- approx(fit_cens$time, fit_cens$surv, xout = Y, rule = 2, yleft = 1)$y
  } else {
    cm <- as.matrix(covars)
    if (!is.null(trt)) cm <- cbind(cm, trt = trt)
    colnames(cm) <- make.names(colnames(cm), unique = TRUE)
    cd <- data.frame(time = time, event = cens_ind, cm)
    tryCatch({
      cf <- as.formula(paste("Surv(time, event) ~", paste(colnames(cm), collapse = "+")))
      cx <- coxph(cf, data = cd, ties = "breslow")
      bh <- basehaz(cx, centered = FALSE)
      lp <- predict(cx, newdata = cd, type = "lp")
      G_hat <- sapply(seq_along(Y), function(i) {
        idx <- max(which(bh$time <= Y[i]))
        if (is.infinite(idx) || idx < 1) 1 else exp(-bh$hazard[idx] * exp(lp[i]))
      })
    }, error = function(e) {
      fit_cens <<- survfit(Surv(time, cens_ind) ~ 1)
      G_hat <<- approx(fit_cens$time, fit_cens$surv, xout = Y, rule = 2, yleft = 1)$y
    })
  }
  G_hat <- pmax(G_hat, eps)
  list(weights = delta_tilde / G_hat, G_hat = G_hat)
}

# =========================================================================
# 2. Per-stratum sorted data
# =========================================================================
make_sorted_stratum <- function(time, status, ipcw, trt, tau) {
  n <- length(time)
  ord <- order(time)
  list(time = time[ord], status = status[ord], ipcw = ipcw[ord], trt = trt[ord],
       ord = ord, iord = order(ord), n = n, tau = tau)
}

# =========================================================================
# 3. Hazard increment influence (same as rmst_cavboost_clean.R)
# =========================================================================
hazard_inc <- function(in_e, in_r, w_total, ipcw_vec, n) {
  denom <- sum(w_total[in_r])
  if (denom <= 0) return(list(h = 0, g = rep(0, n)))
  nume <- sum(w_total[in_e])
  inc <- nume / denom
  g <- (ipcw_vec / denom) * (as.numeric(in_e) - inc * as.numeric(in_r))
  list(h = inc, g = g)
}

# =========================================================================
# 4. Compute per-stratum d_k^(1), d_k^(2) and their gradients
#    All arrays are in the stratum's sorted order
# =========================================================================
stratum_gradient_components <- function(p_sorted, sd, tau) {
  n_k <- sd$n
  q_sorted <- 1 - p_sorted
  
  # Subgroup 1 (p-weighted) and Subgroup 2 ((1-p)-weighted)
  w1 <- p_sorted * sd$ipcw
  w2 <- q_sorted * sd$ipcw
  
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
    ht <- hazard_inc(in_e1, in_r1, w1, sd$ipcw, n_k)
    hc <- hazard_inc(in_e0, in_r0, w1, sd$ipcw, n_k)
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
    ht <- hazard_inc(in_e1, in_r1, w2, sd$ipcw, n_k)
    hc <- hazard_inc(in_e0, in_r0, w2, sd$ipcw, n_k)
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
# 5. Custom XGBoost objective
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
    
    # Get per-stratum d1, d2 and influence gradients g1, g2
    res <- stratum_gradient_components(p_sorted, sd, sd$tau)
    d1 <- res$d1; g1 <- res$g1  # ∂d1/∂p in sorted order
    d2 <- res$d2; g2 <- res$g2  # ∂d2/∂p in sorted order
    
    W_k <- sum(p_k)         # Σ p_i
    n_k_minus_W <- n_k - W_k  # Σ (1-p_i)
    
    # Gradient mirrors original CAVBoost:
    # gp = -(d1 + W·g1 + d2 - (n-W)·g2)
    # grad = p(1-p) · gp  (mapped back to original order)
    # g1 = ∂d1/∂p (p-weighted), g2 = ∂d2/∂(1-p) ((1-p)-weighted)
    # g1/g2 are in SORTED order. sd$iord maps sorted -> original.
    # Extracting g1[sd$iord[j]] gives ∂d1/∂p for original stratum patient j.
    gp_k <- -(d1 + d2 + W_k * g1[sd$iord] + n_k_minus_W * g2[sd$iord])
    grad[idx_k] <- pg[idx_k] * gp_k
  }
  
  list(grad = grad, hess = rep(1.0, n))
}

# =========================================================================
# 6. Prediction
# =========================================================================
pred_stratified <- function(model, dat) {
  features <- setdiff(names(dat), c("trt01p","time","status","A","U","delta","delta_tilde",
                                      "oracle_delta","id","T","C"))
  X <- as.matrix(dat[, features, drop = FALSE])
  1 / (1 + exp(-predict(model, X)))
}

# =========================================================================
# 7. Training
# =========================================================================
train_stratified_cavboost <- function(dat, time, status, tau, stratum,
                                       eta = 0.05, max_depth = 4, nr = 50,
                                       covars = NULL) {
  if ("A" %in% names(dat) && !("trt01p" %in% names(dat))) dat$trt01p <- dat$A
  if ("U" %in% names(dat) && !("time" %in% names(dat))) time <- dat$U
  if ("delta_tilde" %in% names(dat) && !("status" %in% names(dat))) status <- dat$delta_tilde
  
  # IPCW weights (global)
  ipcw_w <- compute_ipcw_weights(time, status, tau, covars, dat$trt01p)
  
  # Pre-sort data per stratum
  K <- max(stratum)
  sorted_strata <- vector("list", K)
  for (k in 1:K) {
    idx <- which(stratum == k)
    if (length(idx) < 2) { sorted_strata[[k]] <- NULL; next }
    sorted_strata[[k]] <- make_sorted_stratum(
      time[idx], status[idx], ipcw_w$weights[idx], dat$trt01p[idx], tau)
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
