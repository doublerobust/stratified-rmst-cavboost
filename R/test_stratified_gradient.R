###############################################################################
# test_stratified_gradient.R
# Verify the per-stratum weighted KM gradient against numerical differentiation.
###############################################################################
set.seed(42)
library(survival); library(xgboost); library(mvtnorm)
source("simulate_data.R"); source("stratified_cavboost.R")
tau <- 12

params <- list(alpha_1=1,alpha_2=1,alpha_3=0.3,alpha_4=0.3,z3_prog=0.3,
               gamma=0,gamma_z3=0,dropout_type="random",
               lambda_drop=NULL,target_censoring_rate=0.30)
dat <- simulate_one_dataset(n=100, params=params)
prog_score <- predict(coxph(Surv(U,delta_tilde)~Z1+Z2+Z3+Z4+Z5+Z6,data=dat),type="lp")
K <- 4
stratum <- as.numeric(cut(prog_score, quantile(prog_score, seq(0,1,1/K)), include.lowest=TRUE))

# Build DMatrix with sorted strata (same as train_stratified_cavboost)
sorted_strata <- vector("list", K)
for (k in 1:K) {
  idx <- which(stratum == k)
  if (length(idx) < 2) next
  sorted_strata[[k]] <- make_sorted_stratum(
    dat$U[idx], dat$delta_tilde[idx], dat$A[idx], tau)
}

Xmat <- as.matrix(dat[, paste0("Z", 1:6)])
dtrain <- xgb.DMatrix(Xmat, label = dat$A)
attr(dtrain, "sorted_strata") <- sorted_strata
attr(dtrain, "stratum") <- stratum

# Analytic gradient at preds=0
set.seed(123)
preds <- rnorm(nrow(dat), 0, 0.5)
analytical <- stratified_loss(preds, dtrain)

# V_strat for numerical differentiation
V_strat <- function(eta) {
  p_raw <- 1/(1+exp(-eta))
  val <- 0
  for (k in 1:K) {
    sd <- sorted_strata[[k]]
    if (is.null(sd) || sd$n < 2) next
    idx_k <- which(stratum == k)
    p_sorted <- p_raw[idx_k][sd$ord]
    res <- stratum_gradient_components(p_sorted, sd, tau)
    W_k <- sum(p_raw[idx_k])
    nk <- length(idx_k)
    val <- val + W_k * res$d1 - (nk - W_k) * res$d2
  }
  val
}

# Numerical gradient
h <- 1e-5
max_diff <- 0
cat("Gradient Verification:\n  j  analytic    numerical   diff\n")
for (j in 1:min(15, nrow(dat))) {
  pu <- pd <- preds; pu[j] <- pu[j] + h; pd[j] <- pd[j] - h
  # numeric = (V(eta+h) - V(eta-h)) / (2h)
  # analytic = -XGBoost_grad (since XGBoost minimizes)
  ana <- -analytical$grad[j]
  num <- (V_strat(pu) - V_strat(pd)) / (2 * h)
  d <- abs(ana - num)
  max_diff <- max(max_diff, d)
  cat(sprintf("  %2d  %10.6f  %10.6f  %10.2e%s\n", j, ana, num, d,
              ifelse(d < 1e-4, " PASS", "")))
}
cat(sprintf("\nMax diff: %.2e  %s\n", max_diff,
            ifelse(max_diff < 1e-4, "PASS", "FAIL")))
