# Run ONE rep of the main comparison (Orig vs Strat_CF vs VT)
# Each rep is a fresh R process, so memory is freed between reps.
# Usage: Rscript run_one_rep_v2.R <scenario> <rep> <outfile>

args <- commandArgs(trailingOnly = TRUE)
sc <- as.numeric(args[1])
rep <- as.numeric(args[2])
outfile <- args[3]

library(mvtnorm); library(ranger)
source("R/stratified_cavboost.R")
source("R/rmst_cavboost_clean.R")

# Manual AUC (avoids pROC/dplyr conflict)
auc_ <- function(p,l) {
  npos <- sum(l); nneg <- sum(!l)
  if (npos < 1 || nneg < 1) return(NA)
  r <- rank(p)
  (sum(r[l]) - npos * (npos + 1) / 2) / (npos * nneg)
}

# DGP
tau <- 30; n_tr <- 500; n_te <- 2000
Smat <- matrix(1/3, 52, 52); diag(Smat) <- 1
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
if (length(unique(oracle)) < 2) { saveRDS(NULL, outfile); q() }

tr <- data.frame(X[1:n_tr, ], trt01p = A[1:n_tr], time = U[1:n_tr], status = st[1:n_tr])
te_df <- data.frame(X[-(1:n_tr), ], trt01p = A[-(1:n_tr)], time = U[-(1:n_tr)], status = st[-(1:n_tr)])
l <- oracle[-(1:n_tr)]

# Original
fo <- tryCatch(train_rmst_cavboost(tr, tr$time, tr$status, tau, eta = 0.05, max_depth = 3, nr = 50), error = function(e) NULL)
po <- if (!is.null(fo)) pred_subgroup(fo, te_df) else rep(0.5, nrow(te_df))
rm(fo); gc()

# Stratified CF (elastic net, all 52 covariates)
feat <- colnames(X)
st_ <- tryCatch(crossfit_prognostic_strata(tr, feat, nfold = 5, K = 4, seed = rep * 100 + sc), error = function(e) NULL)
fs <- if (!is.null(st_)) tryCatch(train_stratified_cavboost(tr, tr$time, tr$status, tau, stratum = st_, eta = 0.1, max_depth = 2, nr = 50), error = function(e) NULL) else NULL
ps <- if (!is.null(fs)) pred_stratified(fs, te_df) else rep(0.5, nrow(te_df))
rm(fs, st_); gc()

# VT (ranger per arm)
ctrl <- tr[tr$trt01p == 0, ]; trt_d <- tr[tr$trt01p == 1, ]
vt <- rep(0, nrow(te_df))
rf_c <- tryCatch(ranger::ranger(Surv(time, status) ~ ., data = ctrl[, !names(ctrl) %in% "trt01p"],
                                 num.trees = 200, min.node.size = 10), error = function(e) NULL)
rf_t <- tryCatch(ranger::ranger(Surv(time, status) ~ ., data = trt_d[, !names(trt_d) %in% "trt01p"],
                                 num.trees = 200, min.node.size = 10), error = function(e) NULL)
if (!is.null(rf_c) && !is.null(rf_t)) {
  pc <- predict(rf_c, te_df[, !names(te_df) %in% "trt01p"])
  pt <- predict(rf_t, te_df[, !names(te_df) %in% "trt01p"])
  tg <- seq(0, tau, length.out = 200)
  surv_grid <- function(s, t, g) { apply(s, 1, function(r) approx(t, r, g, rule = 2, yleft = 1)$y) }
  sc_grid <- surv_grid(pc$survival, pc$unique.death.times, tg)
  st_grid <- surv_grid(pt$survival, pt$unique.death.times, tg)
  vt <- colMeans(st_grid) * tau - colMeans(sc_grid) * tau
}
rm(rf_c, rf_t); gc()

saveRDS(data.frame(sc = sc, rep = rep, orig_auc = auc_(po, l), strat_auc = auc_(ps, l), vt_auc = auc_(vt, l)), outfile)
cat(sprintf("rep %d done\n", rep))
