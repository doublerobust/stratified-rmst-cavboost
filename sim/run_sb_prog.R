set.seed(42)
library(survival); library(xgboost); library(mvtnorm); library(pROC)
source("stratified_cavboost.R"); source("rmst_cavboost_clean.R")
tau <- 12; n_sim <- 50; n <- 600

tg <- seq(0,tau,len=500); dt_tg <- diff(tg)[1]; sh <- 1/0.91
rmst <- function(S) colSums((S[-1,,drop=F]+S[-nrow(S),,drop=F])/2)*dt_tg

run_scenario <- function(prog_effect, label) {
  res <- data.frame()
  for (rep in 1:n_sim) tryCatch({
    # 1 predictive (s1), 3 prognostic (Z1-Z3), 20 noise (z1-z20)
    p_total <- 24
    # All independent for simplicity
    X <- matrix(rnorm((n*2)*p_total), n*2, p_total)
    colnames(X) <- c("s1", paste0("Z",1:3), paste0("z",1:20))
    
    s1 <- X[,"s1"]
    A <- rbinom(n*2, 1, 0.5)
    
    # Treatment effect: moderate, interacts with s1
    # Subgroup: s1 > 0 → benefit, s1 < 0 → no benefit/harm
    trt_effect <- 0.5 + 0.5 * s1
    
    # Prognostic effect (varies by scenario)
    lp_prog <- prog_effect * (X[,"Z1"] + X[,"Z2"] + 0.5*X[,"Z3"])
    
    lp <- 3 + lp_prog + trt_effect * A
    T_lat <- exp(lp + 0.91 * log(-log(runif(n*2))))
    C <- pmin(36 - runif(n*2,0,12), rexp(n*2,rate=0.01))
    U <- pmin(T_lat, C, tau)
    st <- as.numeric(pmin(T_lat, tau) <= C)
    
    # Oracle delta
    S0 <- exp(-outer(tg,exp(3 + lp_prog), function(t,s)(t/s)^sh))
    S1 <- exp(-outer(tg,exp(3 + lp_prog + 0.5 + 0.5*s1), function(t,s)(t/s)^sh))
    od <- rmst(S1) - rmst(S0)
    tb <- as.numeric(od > 0)
    if(length(unique(tb))<2) next
    
    tr_idx <- 1:n; te_idx <- (n+1):(n*2)
    
    dt <- data.frame(X, trt01p=A, time=U, status=st)
    tr <- dt[tr_idx,]; te <- dt[te_idx,]
    tb_tr <- tb[tr_idx]; tb_te <- tb[te_idx]
    
    features <- colnames(X)
    
    # Original CAVBoost
    fit_o <- tryCatch(train_rmst_cavboost(tr,tr$time,tr$status,tau,eta=0.05,max_depth=4,nr=50) ,error=function(e)NULL)
    if(is.null(fit_o)) next
    
    # Stratified CAVBoost (K=4, cross-fitted)
    strat <- tryCatch(crossfit_prognostic_strata(tr, features, nfold=5, K=4, seed=rep*100),
                       error=function(e) NULL)
    fit_s <- if(!is.null(strat)) tryCatch(train_stratified_cavboost(tr,tr$time,tr$status,tau,stratum=strat,eta=0.1,max_depth=3,nr=100) ,error=function(e)NULL) else NULL
    p_s_te <- if(!is.null(fit_s)) pred_stratified(fit_s, te) else rep(0.5,nrow(te))
    
    # VT: Cox per arm
    ctrl <- tr[tr$trt01p==0,]; trt <- tr[tr$trt01p==1,]
    fmla <- paste("Surv(time,status)~", paste(features, collapse="+"))
    fit_c <- suppressWarnings(tryCatch(coxph(as.formula(fmla), data=ctrl, x=TRUE),error=function(e)NULL))
    fit_t <- suppressWarnings(tryCatch(coxph(as.formula(fmla), data=trt, x=TRUE),error=function(e)NULL))
    vt_ite <- if(!is.null(fit_c)&&!is.null(fit_t)) predict(fit_t,newdata=te,type="lp") - predict(fit_c,newdata=te,type="lp") else rep(0,nrow(te))
    
    p_o_te <- pred_subgroup(fit_o, te)
    
    cm <- function(p, lab) {
      r <- suppressMessages(roc(lab, p, quiet=T))
      c(AUC=as.numeric(auc(r)), Acc=mean((p>.5)==lab),
        FPR=if(sum(!lab)>0)sum(p>.5&!lab)/sum(!lab)else NA,
        FNR=if(sum(lab)>0)sum(p<=.5&lab)/sum(lab)else NA)
    }
    
    res <- rbind(res, data.frame(scenario=label,rep=rep,method="OrigCB",t(cm(p_o_te,tb_te))))
    if(!is.null(fit_s)) res <- rbind(res, data.frame(scenario=label,rep=rep,method="StrCB",t(cm(p_s_te,tb_te))))
    res <- rbind(res, data.frame(scenario=label,rep=rep,method="VT",AUC=suppressMessages(roc(tb_te,vt_ite,quiet=T)$auc),Acc=NA,FPR=NA,FNR=NA))
    
    if(rep%%10==0) cat(sprintf("%s rep %d\n",label,rep))
  }, error=function(e) cat(sprintf("skip %s %d: %s\n",label,rep,conditionMessage(e))))
  res
}

# No prognostic effect (pure SubgroupBoost-like)
r1 <- run_scenario(prog_effect=0, "NoProg")

# Strong prognostic effect (Z1=1, Z2=1, Z3=0.5)
r2 <- run_scenario(prog_effect=1, "WithProg")

r <- rbind(r1, r2)
cat("\n========================================\n")
cat("  SUBGROUPBOOST-STYLE DGP (N=600)\n")
cat("========================================\n\n")
for(sc in c("NoProg","WithProg")){
  cat(sprintf("--- %s ---\n",sc))
  cat(sprintf("%-12s %8s %8s %8s %8s %8s\n","Method","AUC","Acc","FPR","FNR","N"))
  for(mn in c("OrigCB","StrCB","VT")){
    s<-r[r$scenario==sc&r$method==mn,]
    cat(sprintf("%-12s %8.4f %8.4f %8.4f %8.4f %8d\n",mn,
        mean(s$AUC,na.rm=T),mean(s$Acc,na.rm=T),mean(s$FPR,na.rm=T),mean(s$FNR,na.rm=T),sum(!is.na(s$AUC))))
  }
}
saveRDS(r,"subgroupboost_prog.rds")
