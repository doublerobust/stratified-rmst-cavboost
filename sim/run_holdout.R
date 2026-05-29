set.seed(42)
library(survival); library(xgboost); library(mvtnorm); library(pROC)
source("simulate_data.R"); source("stratified_cavboost.R"); source("rmst_cavboost_clean.R")
tau <- 12; n_sim <- 50; n_tr <- 300; n_te <- 300

tg <- seq(0,tau,len=500); dt_tg <- diff(tg)[1]; sh <- 1/0.91
rmst <- function(S) colSums((S[-1,,drop=F]+S[-nrow(S),,drop=F])/2)*dt_tg

run_one <- function(params, label) {
  res <- data.frame()
  for (rep in 1:n_sim) tryCatch({
    # Generate
    dat <- simulate_one_dataset(n=n_tr+n_te, params=params)
    if(anyNA(dat$U)) stop("NAs")
    
    # Oracle delta (true subgroup label)
    Zmat <- as.matrix(dat[,paste0("Z",1:6)])
    a0 <- params$alpha_0; a1<-params$alpha_1; a2<-params$alpha_2
    a3<-params$alpha_3; a4<-params$alpha_4; z3p<-params$z3_prog
    S0 <- exp(-outer(tg,exp(a0 + a1*Zmat[,1] + a2*Zmat[,2] + z3p*Zmat[,3]),function(t,s)(t/s)^sh))
    S1 <- exp(-outer(tg,exp(a0 + a1*Zmat[,1] + a2*Zmat[,2] + z3p*Zmat[,3] + a3 + a4*Zmat[,3]),function(t,s)(t/s)^sh))
    od <- rmst(S1)-rmst(S0)
    tb <- as.numeric(od > 0)
    if(length(unique(tb))<2) stop("no var")
    
    # Split
    dat$oracle <- tb
    tr <- dat[1:n_tr,]; te <- dat[(n_tr+1):(n_tr+n_te),]
    
    dt_tr <- tr; names(dt_tr)[names(dt_tr)=="A"]<-"trt01p"
    names(dt_tr)[names(dt_tr)=="U"]<-"time"; names(dt_tr)[names(dt_tr)=="delta_tilde"]<-"status"
    dt_te <- te; names(dt_te)[names(dt_te)=="A"]<-"trt01p"
    names(dt_te)[names(dt_te)=="U"]<-"time"; names(dt_te)[names(dt_te)=="delta_tilde"]<-"status"
    
    # Original CAVBoost
    fit_o <- train_rmst_cavboost(dt_tr,dt_tr$time,dt_tr$status,tau,eta=0.05,max_depth=4,nr=50)
    po_tr <- pred_subgroup(fit_o,dt_tr)
    po_te <- pred_subgroup(fit_o,dt_te)
    
    # Cross-fitted prognostic score → stratified (K=4)
    feat_fit <- c("Z1","Z2","Z3","Z4","Z5","Z6")
    strat <- tryCatch(crossfit_prognostic_strata(dt_tr, feat_fit, nfold=5, K=4, seed=rep*100),
                       error=function(e) NULL)
    fit_s <- if(!is.null(strat)) tryCatch(train_stratified_cavboost(dt_tr,dt_tr$time,dt_tr$status,tau,stratum=strat,eta=0.1,max_depth=3,nr=100) ,error=function(e)NULL) else NULL
    
    if(!is.null(fit_s)){
      ps_tr <- pred_stratified(fit_s,dt_tr)
      ps_te <- pred_stratified(fit_s,dt_te)
    } else { ps_tr <- ps_te <- rep(0.5,n_tr) }
    
    cm <- function(p,t) { r<-roc(t,p,quiet=T)
      list(AUC=as.numeric(auc(r)),Acc=mean((p>.5)==t),
           FPR=if(sum(!t)>0)sum(p>.5&!t)/sum(!t)else NA,
           FNR=if(sum(t)>0)sum(p<=.5&t)/sum(t)else NA) }
    
    for (mn in c("Original","Stratified")) {
      p_tr <- if(mn=="Original") po_tr else ps_tr
      p_te <- if(mn=="Original") po_te else ps_te
      c_tr <- cm(p_tr,tr$oracle); c_te <- cm(p_te,te$oracle)
      res <- rbind(res, data.frame(scenario=label,rep=rep,method=mn,
        AUC_tr=c_tr$AUC,AUC_te=c_te$AUC,Acc_te=c_te$Acc,FPR_te=c_te$FPR,FNR_te=c_te$FNR))
    }
  }, error=function(e) cat(sprintf("skip %s %d\n",label,rep)))
  
  cat(sprintf("%s done (%d reps)\n", label, sum(!is.na(res$AUC_te))))
  res
}

# Base mixed scenario
params_mixed <- list(alpha_0=3,alpha_1=1,alpha_2=1,alpha_3=0.3,alpha_4=0.3,
                     z3_prog=0.3,gamma=0,gamma_z3=0,dropout_type="random",
                     lambda_drop=NULL,target_censoring_rate=0.30)
r1 <- run_one(params_mixed, "Mixed")

# Pure predictive (no prognostic effects)
params_pure <- list(alpha_0=3,alpha_1=0,alpha_2=0,alpha_3=0.3,alpha_4=0.3,
                    z3_prog=0,gamma=0,gamma_z3=0,dropout_type="random",
                    lambda_drop=NULL,target_censoring_rate=0.30)
r2 <- run_one(params_pure, "PurePredictive")

r <- rbind(r1, r2)
cat("\n============================================================\n")
cat("  HOLD-OUT CLASSIFICATION METRICS (N_tr=N_te=300)\n")
cat("============================================================\n\n")
for (sc in c("Mixed","PurePredictive")) {
  cat(sprintf("--- %s ---\n", sc))
  cat(sprintf("%-15s %8s %8s %8s %8s %8s\n","Method","AUC_tr","AUC_te","Acc_te","FPR_te","FNR_te"))
  for(mn in c("Original","Stratified")){
    s<-r[r$scenario==sc&r$method==mn,]
    cat(sprintf("%-15s %8.4f %8.4f %8.4f %8.4f %8.4f\n",mn,
        mean(s$AUC_tr,na.rm=T),mean(s$AUC_te,na.rm=T),
        mean(s$Acc_te,na.rm=T),mean(s$FPR_te,na.rm=T),mean(s$FNR_te,na.rm=T)))
  }
}
saveRDS(r,"holdout_results.rds")
cat("\nSaved to holdout_results.rds\n")
