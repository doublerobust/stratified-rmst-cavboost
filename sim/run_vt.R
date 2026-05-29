set.seed(42)
library(survival); library(xgboost); library(mvtnorm); library(pROC)
source("simulate_data.R"); source("stratified_cavboost.R"); source("rmst_cavboost_clean.R")
tau <- 12; n_sim <- 50; n <- 300

tg <- seq(0,tau,len=500); dt_tg <- diff(tg)[1]; sh <- 1/0.91
rmst_S <- function(S) colSums((S[-1,,drop=F]+S[-nrow(S),,drop=F])/2)*dt_tg

params <- list(alpha_0=3,alpha_1=1,alpha_2=1,alpha_3=0.3,alpha_4=0.3,z3_prog=0.3,
               gamma=0,gamma_z3=0,dropout_type="random",lambda_drop=NULL,target_censoring_rate=0.30)

# VT: Cox per arm, predict RMST at tau
vt_rmst <- function(dat_tr, dat_te) {
  ctrl <- dat_tr[dat_tr$trt01p==0,]
  trt <- dat_tr[dat_tr$trt01p==1,]
  
  fit_c <- tryCatch(coxph(Surv(time,status)~Z1+Z2+Z3+Z4+Z5+Z6,data=ctrl,x=TRUE),error=function(e)NULL)
  fit_t <- tryCatch(coxph(Surv(time,status)~Z1+Z2+Z3+Z4+Z5+Z6,data=trt,x=TRUE),error=function(e)NULL)
  if(is.null(fit_c)||is.null(fit_t)) return(rep(NA,nrow(dat_te)))
  
  lp_c <- predict(fit_c, newdata=dat_te, type="lp")
  lp_t <- predict(fit_t, newdata=dat_te, type="lp")
  
  rmst_from_cox <- function(lp, bh, tau) {
    # Build cumulative hazard at event times
    H0 <- bh$hazard
    H0_times <- bh$time
    # For each patient
    rmst <- numeric(length(lp))
    for(i in seq_along(lp)) {
      # H(t) = H0(t) * exp(lp[i]) at each event time
      H <- H0 * exp(lp[i])
      S <- exp(-H)
      # All event times up to tau, plus endpoints
      idx <- which(H0_times <= tau)
      if(length(idx)==0) { rmst[i] <- tau; next }
      t_vals <- c(0, H0_times[idx], tau)
      s_vals <- c(1, S[idx], if(max(idx)<length(S)) S[max(idx)] else 1)
      # Remove duplicate tau if present
      if(t_vals[length(t_vals)-1]==tau) { t_vals <- t_vals[-length(t_vals)]; s_vals <- s_vals[-length(s_vals)] }
      rmst[i] <- sum(diff(t_vals) * s_vals[-length(s_vals)])
    }
    rmst
  }
  
  bh_c <- basehaz(fit_c, centered=FALSE)
  bh_t <- basehaz(fit_t, centered=FALSE)
  rmst_c <- rmst_from_cox(lp_c, bh_c, tau)
  rmst_t <- rmst_from_cox(lp_t, bh_t, tau)
  
  rmst_t - rmst_c
}

results <- data.frame()
for (rep in 1:n_sim) tryCatch({
  dat <- simulate_one_dataset(n=n*2, params=params)
  if(anyNA(dat$U)) stop("NAs")
  
  Z <- as.matrix(dat[,paste0("Z",1:6)])
  S0 <- exp(-outer(tg,exp(3+Z[,1]+Z[,2]+0.3*Z[,3]),function(t,s)(t/s)^sh))
  S1 <- exp(-outer(tg,exp(3+Z[,1]+Z[,2]+0.3*Z[,3]+0.3+0.3*Z[,3]),function(t,s)(t/s)^sh))
  od <- rmst_S(S1)-rmst_S(S0)
  tb <- as.numeric(od > 0)
  
  tr <- 1:n; te <- (n+1):(n*2)
  dt_tr <- dat[tr,]; names(dt_tr)[names(dt_tr)=="A"]<-"trt01p"
  names(dt_tr)[names(dt_tr)=="U"]<-"time"; names(dt_tr)[names(dt_tr)=="delta_tilde"]<-"status"
  dt_te <- dat[te,]; names(dt_te)[names(dt_te)=="A"]<-"trt01p"
  names(dt_te)[names(dt_te)=="U"]<-"time"; names(dt_te)[names(dt_te)=="delta_tilde"]<-"status"
  
  vt <- vt_rmst(dt_tr, dt_te)
  
  fit_o <- train_rmst_cavboost(dt_tr,dt_tr$time,dt_tr$status,tau,eta=0.05,max_depth=4,nr=50)
  po_te <- pred_subgroup(fit_o,dt_te)
  
  strat <- tryCatch(crossfit_prognostic_strata(dt_tr, paste0("Z",1:6), nfold=5, K=4, seed=label*100),
                     error=function(e) NULL)
  fit_s <- if(!is.null(strat)) tryCatch(train_stratified_cavboost(dt_tr,dt_tr$time,dt_tr$status,tau,stratum=strat,eta=0.1,max_depth=3,nr=100) ,error=function(e)NULL) else NULL
  ps_te <- if(!is.null(fit_s)) pred_stratified(fit_s,dt_te) else rep(0.5,n)
  
  cm <- function(p,t) {
    r<-roc(t,p,quiet=T)
    c(AUC=as.numeric(auc(r)),Acc=mean((p>.5)==t),FPR=if(sum(!t)>0)sum(p>.5&!t)/sum(!t)else NA,FNR=if(sum(t)>0)sum(p<=.5&t)/sum(t)else NA)
  }
  
  results <- rbind(results, data.frame(rep=rep,method="VT",t(c(cm(vt,tb[te]),sd_p=sd(vt,na.rm=T)))))
  results <- rbind(results, data.frame(rep=rep,method="OrigCB",t(c(cm(po_te,tb[te]),sd_p=sd(po_te)))))
  if(!is.null(fit_s)) results <- rbind(results, data.frame(rep=rep,method="StrCB",t(c(cm(ps_te,tb[te]),sd_p=sd(ps_te)))))
}, error=function(e) cat(sprintf("skip %d: %s\n",rep,conditionMessage(e))))

cat("\n=== VT (Cox RMST) vs CAVBoost (hold-out N=300) ===\n")
cat(sprintf("%-10s %8s %8s %8s %8s %8s %8s\n","Method","AUC","Acc","FPR","FNR","sd","N"))
for(mn in c("VT","OrigCB","StrCB")){
  s<-results[results$method==mn,]
  cat(sprintf("%-10s %8.4f %8.4f %8.4f %8.4f %8.4f %8d\n",mn,
      mean(s$AUC,na.rm=T),mean(s$Acc,na.rm=T),mean(s$FPR,na.rm=T),mean(s$FNR,na.rm=T),
      mean(s$sd_p,na.rm=T),sum(!is.na(s$AUC))))
}
saveRDS(results,"vt_comparison.rds")
