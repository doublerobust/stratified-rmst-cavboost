set.seed(42)
library(survival); library(xgboost); library(mvtnorm); library(pROC); library(pseudo)
source("simulate_data.R"); source("stratified_cavboost.R"); source("rmst_cavboost_clean.R")
tau <- 12; n_sim <- 50
params <- list(alpha_1=1,alpha_2=1,alpha_3=0.3,alpha_4=0.3,z3_prog=0.3,gamma=0,gamma_z3=0,dropout_type="random",lambda_drop=NULL,target_censoring_rate=0.30)

cm <- function(p,t) {
  r <- tryCatch(roc(t,p,quiet=T),error=function(e)NULL)
  if(is.null(r)) return(c(AUC=NA,Acc=NA,FPR=NA,FNR=NA))
  c(AUC=as.numeric(auc(r)),Acc=mean((p>.5)==t),FPR=if(sum(!t)>0)sum(p>.5&!t)/sum(!t)else NA,FNR=if(sum(t)>0)sum(p<=.5&t)/sum(t)else NA)
}

prog_score <- function(dt, method) {
  switch(method,
    CoxPH = predict(coxph(Surv(time,status)~Z1+Z2+Z3+Z4+Z5+Z6,data=dt),type="lp"),
    PseudoXGB = {
      pr <- pseudomean(dt$time,dt$status,tmax=tau)
      X <- as.matrix(dt[,paste0("Z",1:6)])
      predict(xgb.train(params=list(objective="reg:squarederror",eta=0.1,max_depth=2),data=xgb.DMatrix(X,label=pr),nrounds=50,verbose=0),X)
    },
    XGBCox = {
      X <- as.matrix(dt[,paste0("Z",1:6)]); dtm<-xgb.DMatrix(X)
      setinfo(dtm,"label",dt$time); setinfo(dtm,"weight",dt$status)
      predict(xgb.train(params=list(objective="survival:cox",eta=0.1,max_depth=2),data=dtm,nrounds=50,verbose=0),X)
    }
  )
}

K_vals <- c(2,4)
prog_names <- c("CoxPH","PseudoXGB","XGBCox")
res <- data.frame()
for (rep in 1:n_sim) tryCatch({
  dat <- simulate_one_dataset(n=300,params=params)
  if(anyNA(dat$U)) stop("NAs")
  tb <- as.numeric(dat$oracle_delta>0)
  if(length(unique(tb))<2) stop("no var")
  dt<-dat;names(dt)[names(dt)=="A"]<-"trt01p"
  names(dt)[names(dt)=="U"]<-"time";names(dt)[names(dt)=="delta_tilde"]<-"status"
  
  fit_o <- train_rmst_cavboost(dt,dt$time,dt$status,tau,eta=0.05,max_depth=4,nr=50,covars=dt[,paste0("Z",1:6)])
  p_o <- pred_subgroup(fit_o,dt)
  res <- rbind(res, data.frame(rep=rep,method="Original",prog="none",K=NA,t(cm(p_o,tb))))
  
  for(pn in prog_names){
    ps <- tryCatch(prog_score(dt,pn),error=function(e)NULL)
    if(is.null(ps)) next
    for(Ki in K_vals){
      strat <- as.numeric(cut(ps,quantile(ps,seq(0,1,1/Ki),na.rm=T),include.lowest=TRUE))
      fit_s <- tryCatch(train_stratified_cavboost(dt,dt$time,dt$status,tau,stratum=strat,eta=0.1,max_depth=3,nr=100,covars=dt[,paste0("Z",1:6)]),error=function(e)NULL)
      if(!is.null(fit_s)) res <- rbind(res, data.frame(rep=rep,method="Stratified",prog=pn,K=Ki,t(cm(pred_stratified(fit_s,dt),tb))))
    }
  }
  if(rep%%10==0) cat(sprintf("rep %d/%d\n",rep,n_sim))
}, error=function(e) cat(sprintf("skip %d\n",rep)))

cat("\n=== Classification metrics by prognostic score (N=300) ===\n")
for(mn in c("Original","Stratified")){
  if(mn=="Original"){
    s<-res[res$method=="Original",]
    cat(sprintf("\nOriginal (no stratification): AUC=%.4f Acc=%.4f FPR=%.4f FNR=%.4f\n",
                mean(s$AUC,na.rm=T),mean(s$Acc,na.rm=T),mean(s$FPR,na.rm=T),mean(s$FNR,na.rm=T)))
  } else {
    cat("\nStratified:\n")
    cat(sprintf("%-12s %4s %8s %8s %8s %8s\n","ProgScore","K","AUC","Acc","FPR","FNR"))
    s<-res[res$method=="Stratified",]
    for(pn in prog_names) for(Ki in K_vals){
      ss<-s[s$prog==pn&s$K==Ki,]
      cat(sprintf("%-12s  %2d %8.4f %8.4f %8.4f %8.4f\n",pn,Ki,mean(ss$AUC,na.rm=T),mean(ss$Acc,na.rm=T),mean(ss$FPR,na.rm=T),mean(ss$FNR,na.rm=T)))
    }
  }
}
