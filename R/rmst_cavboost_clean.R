###############################################################################
# rmst_cavboost_clean.R
# RMST-CAVBoost with weighted Kaplan-Meier and ANALYTIC gradients.
#
# Architecture (same as SubgroupBoost.RMST):
#   Weighted KM per arm x subgroup -> RMST -> dRMST -> Loss -> Analytic Gradient
#
# Key improvements over rmst_cavboost_km.R:
#   1. Sort data by event time ONCE before the boosting loop
#   2. inc_fn defined OUTSIDE the for loop (avoids scope issues)
#   3. gH accumulates as VECTOR of length n (not per-sg index)
#   4. gRMST accumulates as VECTOR of length n
#   5. Explicit sign flip for sg2 gradients (weights = 1-p)
#   6. All operations vectorized (no sapply, no numerical gradients)
###############################################################################

library(xgboost)
library(survival)
library(dplyr)


make_sorted_data <- function(time, status, trt, tau) {
  n <- length(time)
  ord <- order(time)
  list(time = time[ord], status = status[ord], trt = trt[ord],
       ord = ord, iord = order(ord), n = n, tau = tau)
}

hazard_inc <- function(in_e, in_r, w_total, n) {
  denom <- sum(w_total[in_r])
  if (denom <= 0) return(list(h = 0, g = rep(0, n)))
  nume <- sum(w_total[in_e])
  inc <- nume / denom
  g <- (1 / denom) * (as.numeric(in_e) - inc * as.numeric(in_r))
  list(h = inc, g = g)
}

rmst_cavboost_loss <- function(preds, dtrain) {
  trt <- getinfo(dtrain, "label")
  sd <- attr(dtrain, "sorted_data")
  tau <- sd$tau
  n <- length(preds)
  p_raw <- 1 / (1 + exp(-preds))
  pg <- exp(preds) / (1 + exp(preds))^2
  p <- p_raw[sd$ord]
  event_times <- unique(sd$time[sd$status == 1])
  event_times <- event_times[event_times <= tau]
  if (length(event_times) == 0) return(list(grad = rep(0, n), hess = rep(0.0001, n)))
  dt <- event_times - c(0, event_times[-length(event_times)])
  av <- c(1, 0)
  H <- c(0, 0, 0, 0)
  gH <- list(rep(0, n), rep(0, n), rep(0, n), rep(0, n))
  dRMST <- c(0, 0)
  gRMST <- list(rep(0, n), rep(0, n))
  w1 <- p
  w2 <- (1 - p)
  for (idx in seq_along(event_times)) {
    ti <- event_times[idx]
    dti <- dt[idx]
    in_r1 <- sd$time >= ti & sd$trt == av[1]
    in_r0 <- sd$time >= ti & sd$trt == av[2]
    in_e1 <- sd$time == ti & sd$status == 1 & sd$trt == av[1]
    in_e0 <- sd$time == ti & sd$status == 1 & sd$trt == av[2]
    r1t <- hazard_inc(in_e1, in_r1, w1, n)
    r1c <- hazard_inc(in_e0, in_r0, w1, n)
    r2t <- hazard_inc(in_e1, in_r1, w2, n)
    r2c <- hazard_inc(in_e0, in_r0, w2, n)
    H <- H + c(r1t$h, r1c$h, r2t$h, r2c$h)
    gH[[1]] <- gH[[1]] + r1t$g
    gH[[2]] <- gH[[2]] + r1c$g
    gH[[3]] <- gH[[3]] - r2t$g
    gH[[4]] <- gH[[4]] - r2c$g
    S <- exp(-H)
    dRMST[1] <- dRMST[1] + (S[1] - S[2]) * dti
    dRMST[2] <- dRMST[2] + (S[3] - S[4]) * dti
    gRMST[[1]] <- gRMST[[1]] + (-S[1] * gH[[1]] + S[2] * gH[[2]]) * dti
    gRMST[[2]] <- gRMST[[2]] + (-S[3] * gH[[3]] + S[4] * gH[[4]]) * dti
  }
  sp <- sum(p_raw)
  sq <- sum(1 - p_raw)
  gp <- -(dRMST[1] + sp * gRMST[[1]] + dRMST[2] - sq * gRMST[[2]])
  grad <- as.numeric(pg * gp[sd$iord])
  list(grad = grad, hess = rep(1.0, n))
}

rmst_cavboost_eval <- function(preds, dtrain) {
  trt <- getinfo(dtrain, "label")
  sd <- attr(dtrain, "sorted_data")
  tau <- sd$tau
  p_raw <- 1 / (1 + exp(-preds))
  p <- p_raw[sd$ord]
  av <- c(1, 0)
  event_times <- unique(sd$time[sd$status == 1])
  event_times <- event_times[event_times <= tau]
  if (length(event_times) == 0) return(list(metric = "OTR_error", value = 0))
  dt <- event_times - c(0, event_times[-length(event_times)])
  H <- c(0, 0, 0, 0)
  dRMST <- c(0, 0)
  w1 <- p
  w2 <- (1 - p)
  for (idx in seq_along(event_times)) {
    ti <- event_times[idx]; dti <- dt[idx]
    in_r1 <- sd$time >= ti & sd$trt == av[1]
    in_r0 <- sd$time >= ti & sd$trt == av[2]
    in_e1 <- sd$time == ti & sd$status == 1 & sd$trt == av[1]
    in_e0 <- sd$time == ti & sd$status == 1 & sd$trt == av[2]
    h1 <- function(in_e, in_r, w) { d <- sum(w[in_r]); if (d <= 0) 0 else sum(w[in_e]) / d }
    H <- H + c(h1(in_e1, in_r1, w1), h1(in_e0, in_r0, w1),
               h1(in_e1, in_r1, w2), h1(in_e0, in_r0, w2))
    S <- exp(-H)
    dRMST[1] <- dRMST[1] + (S[1] - S[2]) * dti
    dRMST[2] <- dRMST[2] + (S[3] - S[4]) * dti
  }
  err <- -(sum(p_raw) * dRMST[1] - sum(1 - p_raw) * dRMST[2])
  list(metric = "OTR_error", value = if (is.na(err) || is.infinite(err)) 1e10 else err)
}
train_rmst_cavboost <- function(dat, time, status, tau,
                                 eta = 0.05, max_depth = 4, nr = 50, covars = NULL) {
  sorted_data <- make_sorted_data(time, status, dat$trt01p, tau)
  Xmat <- as.matrix(select(dat, -trt01p, -time, -status))
  dtrain <- xgb.DMatrix(Xmat, label = dat$trt01p)
  attr(dtrain, "sorted_data") <- sorted_data
  xgb.train(params = list(eta = eta, max_depth = max_depth, lambda = 1,
                           min_child_weight = 0, subsample = 1, colsample_bytree = 1),
            data = dtrain, nrounds = nr,
            objective = rmst_cavboost_loss, custom_metric = rmst_cavboost_eval, verbose = 0)
}

pred_subgroup <- function(model, dat) {
  1 / (1 + exp(-predict(model, as.matrix(select(dat, -trt01p, -time, -status)))))
}

gen_surv_data <- function(N) {
  Z1 <- rnorm(N); Z2 <- rnorm(N); Z3 <- rnorm(N); Z4 <- rnorm(N)
  Z5 <- rnorm(N); Z6 <- rnorm(N); Z7 <- rnorm(N)
  A <- rbinom(N, 1, 0.5)
  sg <- as.numeric(Z7 > 0)
  logHR <- ifelse(A == 1, ifelse(sg == 1, log(0.2), log(3.0)), 0)
  eta <- 0.3 * Z1 + 0.3 * Z2 + 0.3 * Z3 + logHR
  shape <- 1.2; scale <- 8.0
  T_latent <- (-log(runif(N)) / (1/scale^shape * exp(eta)))^(1/shape)
  tau <- quantile(T_latent, 0.7)
  C <- runif(N, tau * 0.5, tau * 1.5)
  time <- pmin(T_latent, C)
  status <- as.numeric(T_latent <= C)
  list(dat = data.frame(Z1, Z2, Z3, Z4, Z5, Z6, Z7, trt01p = A, time = time, status = status),
       true = sg, tau = as.numeric(tau))
}

# Test suite (always runs when sourced directly)
if (interactive() || Sys.getenv("RMST_CLEAN_TEST") == "1") {
  library(pROC)
  set.seed(42)
  cat("=== RMST-CAVBoost (Clean Implementation) ===\n\n")
  cat("Generating data...\n")
  train <- gen_surv_data(600)
  test <- gen_surv_data(2000)
  cat(sprintf("Train: n=%d, true sg proportion: %.3f\n", nrow(train$dat), mean(train$true)))
  cat("\n--- Gradient Verification ---\n")
  cat("Comparing analytic vs numerical gradient (first 5 obs)...\n")
  sub_idx <- 1:200
  mtest <- train$dat[sub_idx, ]
  sd_sub <- make_sorted_data(mtest$time, mtest$status, mtest$trt01p, train$tau)
  Xsub <- as.matrix(select(mtest, -trt01p, -time, -status))
  dt_sub <- xgb.DMatrix(Xsub, label = mtest$trt01p)
  attr(dt_sub, "sorted_data") <- sd_sub
  preds0 <- rep(0, 200)
  r <- rmst_cavboost_loss(preds0, dt_sub)
  h <- 0.001; max_diff <- 0
  for (j in 1:5) {
    pu <- preds0; pu[j] <- pu[j] + h
    pd <- preds0; pd[j] <- pd[j] - h
    nu <- rmst_cavboost_eval(pu, dt_sub)$value
    nd <- rmst_cavboost_eval(pd, dt_sub)$value
    num_grad <- (nu - nd) / (2 * h)
    diff <- abs(r$grad[j] - num_grad)
    max_diff <- max(max_diff, diff)
    cat(sprintf("  j=%d: analytic=%.6f  numerical=%.6f  diff=%.2e\n", j, r$grad[j], num_grad, diff))
  }
  cat(sprintf("  Max diff: %.2e  %s\n\n", max_diff, ifelse(max_diff < 1e-4, "PASS", "FAIL")))
  cat("--- Training ---\n")
  for (nr in c(5, 10, 20)) {
    cat(sprintf("\nnr = %d rounds:\n", nr))
    tm <- system.time({
      m <- train_rmst_cavboost(dat = train$dat, time = train$dat$time,
        status = train$dat$status, tau = train$tau, eta = 0.05, max_depth = 4, nr = nr,
        covars = train$dat[, c("Z1","Z2","Z3","Z4","Z5","Z6","Z7")])
    })
    pp <- pred_subgroup(m, test$dat)
    pc <- as.numeric(pp > 0.5); t <- test$true
    acc <- mean(pc == t)
    roc_obj <- roc(t, pp, quiet = TRUE); auc_val <- auc(roc_obj)
    imp <- tryCatch(xgb.importance(model = m,
      feature_names = setdiff(colnames(train$dat), c("trt01p","time","status"))), error = function(e) NULL)
    top_var <- if (!is.null(imp) && nrow(imp) > 0) as.character(imp$Feature[1]) else "?"
    cat(sprintf("  Time: %.1fs  Acc=%.3f  AUC=%.3f  p=[%.3f, %.3f]  Top=%s\n",
      tm[3], acc, auc_val, min(pp), max(pp), top_var))
    if (!is.null(imp) && nrow(imp) > 0) {
      cat("  Top 3 vars:\n")
      print(imp[1:min(3, nrow(imp)), ])
    }
  }
  cat("\nDone.\n")
}
