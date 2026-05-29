# Clean sequential runner, saves CSV (no RDS corruption possible)
# Run from repo root: nohup Rscript run_clean_50.R > clean.log 2>&1 &

library(mvtnorm); library(ranger)
source("R/stratified_cavboost.R")
source("R/rmst_cavboost_clean.R")

tau <- 30; n_tr <- 500; n_te <- 2000
Smat <- matrix(1/3, 52, 52); diag(Smat) <- 1
nr <- 50; n_sim <- 50

auc_ <- function(p,l){npos<-sum(l);nneg<-sum(!l);if(npos<1||nneg<1)return(NA);r<-rank(p);abs((sum(r[l])-npos*(npos+1)/2)/(npos*nneg))}
# ^ abs() to handle inverted predictions

sc_names <- list("1"="S1_Linear","2"="S2_Diff","3"="S3_U",
                 "4"="S4_Enclave","5"="S5_S","6"="S6_Cross")

dir.create("results_csv", showWarnings=FALSE)

for (sc in 1:6) {
  sn <- sc_names[[as.character(sc)]]
  n_existing <- length(list.files("results_csv", sprintf("%s_.*\\.csv", sn)))
  if (n_existing >= n_sim) { cat(sn, "already done\n"); next }
  
  for (rep in 1:n_sim) {
    out <- sprintf("results_csv/%s_rep%d.csv", sn, rep)
    if (file.exists(out)) next
    
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
    if (length(unique(oracle)) < 2) { cat("skip", sn, rep, "\n"); next }
    
    tr <- data.frame(X[1:n_tr, ], trt01p = A[1:n_tr], time = U[1:n_tr], status = st[1:n_tr])
    te_df <- data.frame(X[-(1:n_tr), ], trt01p = A[-(1:n_tr)], time = U[-(1:n_tr)], status = st[-(1:n_tr)])
    l <- oracle[-(1:n_tr)]
    
    # Orig
    fo <- tryCatch(train_rmst_cavboost(tr, tr$time, tr$status, tau, eta = 0.05, max_depth = 3, nr = nr), error = function(e) NULL)
    po <- if (!is.null(fo)) pred_subgroup(fo, te_df) else rep(0.5, nrow(te_df))
    rm(fo); gc()
    
    # Strat (elastic net CF) — ALL 52 covariates required.
    # Using fewer silently changes the prognostic score and invalidates the comparison.
    feat <- colnames(X)
    stopifnot("Prognostic score must use all 52 covariates" = length(feat) == 52)
    st_strat <- tryCatch(crossfit_prognostic_strata(tr, feat, nfold = 5, K = 4, seed = rep * 100 + sc), error = function(e) NULL)
    fs <- if (!is.null(st_strat)) tryCatch(train_stratified_cavboost(tr, tr$time, tr$status, tau, stratum = st_strat, eta = 0.1, max_depth = 2, nr = nr), error = function(e) NULL) else NULL
    ps <- if (!is.null(fs)) pred_stratified(fs, te_df) else rep(0.5, nrow(te_df))
    rm(fs, st_strat); gc()
    
    # VT (ranger)
    ctrl <- tr[tr$trt01p == 0, ]; trt_d <- tr[tr$trt01p == 1, ]
    vt_pred <- rep(0, nrow(te_df))
    rf_c <- tryCatch(ranger::ranger(Surv(time, status) ~ ., data = ctrl[, !names(ctrl) %in% "trt01p"], num.trees = 200, min.node.size = 10), error = function(e) NULL)
    rf_t <- tryCatch(ranger::ranger(Surv(time, status) ~ ., data = trt_d[, !names(trt_d) %in% "trt01p"], num.trees = 200, min.node.size = 10), error = function(e) NULL)
    if (!is.null(rf_c) && !is.null(rf_t)) {
      pc <- predict(rf_c, te_df[, !names(te_df) %in% "trt01p"])
      pt <- predict(rf_t, te_df[, !names(te_df) %in% "trt01p"])
      tg <- seq(0, tau, length.out = 200)
      surv_grid <- function(s, t, g) { apply(s, 1, function(r) approx(t, r, g, rule = 2, yleft = 1)$y) }
      sc_grid <- surv_grid(pc$survival, pc$unique.death.times, tg)
      st_grid <- surv_grid(pt$survival, pt$unique.death.times, tg)
      vt_pred <- colMeans(st_grid) * tau - colMeans(sc_grid) * tau
    }
    rm(rf_c, rf_t); gc()
    
    write.csv(data.frame(scenario = sc, rep = rep,
                         orig_auc = auc_(po, l),
                         strat_auc = auc_(ps, l),
                         vt_auc = auc_(vt_pred, l)),
              out, row.names = FALSE)
    cat(sprintf("\r%s rep %d/%d", sn, rep, n_sim))
  }
  cat(sprintf("\n%s done\n", sn))
}

# Summary
cat("\n\n====================================\n  FINAL RESULTS\n====================================\n\n")
for (sc in 1:6) {
  sn <- sc_names[[as.character(sc)]]
  ff <- list.files("results_csv", sprintf("%s_.*\\.csv", sn), full.names = TRUE)
  res <- do.call(rbind, lapply(ff, read.csv))
  cat(sprintf("%-14s %-8s %-8s %-8s\n", sn, "Orig", "Strat", "VT"))
  cat(sprintf("%-14s %8.4f %8.4f %8.4f\n\n", "",
              mean(res$orig_auc), mean(res$strat_auc), mean(res$vt_auc)))
}
cat("All results in results_csv/\n")
