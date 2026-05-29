# Holdout model selection: train both, pick by held-out RMST gain
# Usage: Rscript run_holdout_selection.R (from repo root)

set.seed(42)
library(survival); library(xgboost); library(mvtnorm); library(pROC)
source("R/stratified_cavboost.R")
source("~/.openclaw/workspace/CAVBoost/rmst_cavboost_clean.R")

tau <- 30; n_sim <- 30; n_tr <- 500; n_te <- 2000; rho <- 1/3

gen_data <- function(scenario, prog_setting=2, n=500) {
  p_total <- 52
  Sigma <- matrix(rho, p_total, p_total); diag(Sigma) <- 1
  X <- rmvnorm(n, sigma = Sigma)
  colnames(X) <- c(paste0("z",1:50), "S1", "S2")
  A <- rbinom(n, 1, 0.5); beta0 <- sqrt(6); sigma0 <- 0.4
  beta_Z <- if(prog_setting==1) rep(0,50) else c(rep(0.4,4), rep(0,46))
  Zbeta <- X[,1:50] %*% beta_Z
  treat_eff <- switch(scenario,
    `1` = X[,"S1"],
    `3` = { s<-X[,"S1"]; 2*ifelse(abs(s)<0.67,exp(-s^2)-0.4,exp(-s^2)-0.8) },
    `4` = { s1<-X[,"S1"];s2<-X[,"S2"]; 2*((-1.07<=s1&s1<1.07)&(-1.07<=s2&s2<1.07))-1 },
    `5` = { s<-X[,"S1"]; 2*ifelse(s>=0.67|(-0.67<=s&s<0),1,0)-1 }
  )
  if (scenario %in% 4:6) Zbeta <- -Zbeta^2
  T_lat <- exp(beta0 + A*treat_eff + Zbeta + sigma0*rnorm(n))
  C <- pmin(30, rexp(n, rate=-log(0.9)/12))
  U <- pmin(T_lat, C, tau)
  st <- as.numeric(T_lat <= C)
  list(X=X, A=A, U=U, st=st, oracle=as.numeric(treat_eff>0), tau=tau)
}

run_sc <- function(sc, label) {
  res <- data.frame()
  for (rep in 1:n_sim) tryCatch({
    d <- gen_data(sc, 2, n_tr+n_te)
    if(length(unique(d$oracle))<2) next
    tr <- data.frame(d$X[1:n_tr,], trt01p=d$A[1:n_tr], time=d$U[1:n_tr], status=d$st[1:n_tr])
    te <- data.frame(d$X[-(1:n_tr),], trt01p=d$A[-(1:n_tr)], time=d$U[-(1:n_tr)], status=d$st[-(1:n_tr)])
    l_te <- d$oracle[-(1:n_tr)]
    f <- colnames(d$X)

    # Original
    fo <- tryCatch(train_rmst_cavboost(tr,tr$time,tr$status,tau,eta=0.05,max_depth=3,nr=50),error=function(e)NULL)
    po <- if(!is.null(fo)) pred_subgroup(fo,te) else rep(0.5,nrow(te))

    # Stratified (cross-fitted)
    # Cross-fit on 6 prognostic vars (not all 52 — faster and avoids overfit)
    feat_prog <- c("Z1","Z2","Z3","Z4","S1","S2")
    st_ <- tryCatch(crossfit_prognostic_strata(tr, feat_prog, nfold=5, K=4, seed=rep*100+sc), error=function(e) NULL)
    ps_te <- if(!is.null(st_)) {
      fs <- tryCatch(train_stratified_cavboost(tr,tr$time,tr$status,tau,stratum=st_,eta=0.1,max_depth=2,nr=50),error=function(e)NULL)
      if(!is.null(fs)) pred_stratified(fs,te) else rep(0.5,nrow(te))
    } else rep(0.5,nrow(te))

    # Holdout selection (stratified + original, pick by RMST gain)
    sel <- tryCatch(select_model_by_holdout(tr, "time", "status", tau, stratum=st_, features=feat_prog),
                     error=function(e) NULL)
    if (!is.null(sel)) {
      psel_te <- if (sel$winner == "original") pred_subgroup(sel$fit, te) else pred_stratified(sel$fit, te)
    } else psel_te <- rep(0.5, nrow(te))

    cmf <- function(p,l) {
      r<-suppressMessages(roc(l,p,quiet=T))
      c(AUC=as.numeric(auc(r)),Acc=mean((p>.5)==l),
        FPR=if(sum(!l)>0)sum(p>.5&!l)/sum(!l)else NA,FNR=if(sum(l)>0)sum(p<=.5&l)/sum(l)else NA)
    }
    res <- rbind(res, data.frame(sc=label,rep=rep,method="OrigCB",t(cmf(po,l_te))))
    res <- rbind(res, data.frame(sc=label,rep=rep,method="StrCB_CF",t(cmf(ps_te,l_te))))
    res <- rbind(res, data.frame(sc=label,rep=rep,method="Select",t(cmf(psel_te,l_te))))
    if (!is.null(sel)) {
      attr(res, "winners") <- rbind(attr(res, "winners", exact=TRUE),
        data.frame(sc=label, rep=rep, winner=sel$winner, gain_orig=sel$gain_orig, gain_strat=sel$gain_strat))
    }
  }, error=function(e) cat(sprintf("skip %s %d: %s\n",label,rep,e$message)))
  cat(sprintf("%s done (%d reps)\n",label,sum(!is.na(res$AUC))))
  res
}

all <- data.frame()
for (sc in list(c(1,"S1_Linear"),c(3,"S3_U"),c(4,"S4_Enclave"),c(5,"S5_S"))) {
  all <- rbind(all, run_sc(as.numeric(sc[1]),sc[2]))
}

cat("\n============================================================\n")
cat("  HOLDOUT MODEL SELECTION (30 reps, CF strata)\n")
cat("============================================================\n")
for (sc in c("S1_Linear","S3_U","S4_Enclave","S5_S")) {
  cat(sprintf("\n--- %s ---\n", sc))
  cat(sprintf("%-12s %8s %8s %8s %8s\n","Method","AUC","Acc","FPR","FNR"))
  for(mn in c("OrigCB","StrCB_CF","Select")){
    s<-all[all$sc==sc&all$method==mn,]
    cat(sprintf("%-12s %8.4f %8.4f %8.4f %8.4f\n",mn,
        mean(s$AUC,na.rm=T),mean(s$Acc,na.rm=T),mean(s$FPR,na.rm=T),mean(s$FNR,na.rm=T)))
  }
  # Selection rate
  w<-attr(all, "winners", exact=TRUE)
  if (!is.null(w)) {
    w_sc <- w[w$sc==sc,]
    cat(sprintf("  Selected stratified: %.0f/%.0f reps\n",
                sum(w_sc$winner=="stratified",na.rm=T), nrow(w_sc)))
  }
}

saveRDS(all, "holdout_selection_results.rds")
cat("\nSaved to holdout_selection_results.rds\n")
