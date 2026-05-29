set.seed(42)
library(survival); library(xgboost); library(mvtnorm); library(pROC)
source("stratified_cavboost.R"); source("rmst_cavboost_clean.R")
tau <- 30; n_sim <- 30; n_tr <- 500; n_te <- 2000; rho <- 1/3

gen_data <- function(scenario, prog_setting=2, n=500) {
  p_total <- 52
  Sigma <- matrix(rho, p_total, p_total); diag(Sigma) <- 1
  X <- rmvnorm(n, sigma = Sigma)
  colnames(X) <- c(paste0("z",1:50), "S1", "S2")
  A <- rbinom(n, 1, 0.5)
  beta0 <- sqrt(6); sigma0 <- 0.4
  beta_Z <- if(prog_setting==1) rep(0,50) else c(rep(0.4,4), rep(0,46))
  Zbeta <- X[,1:50] %*% beta_Z

  treat_eff <- switch(scenario,
    `1` = X[,"S1"],
    `2` = X[,"S1"] - X[,"S2"],
    `3` = { s<-X[,"S1"]; 2*ifelse(abs(s)<0.67,exp(-s^2)-0.4,exp(-s^2)-0.8) },
    `4` = { s1<-X[,"S1"];s2<-X[,"S2"]; 2*((-1.07<=s1&s1<1.07)&(-1.07<=s2&s2<1.07))-1 },
    `5` = { s<-X[,"S1"]; 2*ifelse(s>=0.67|(-0.67<=s&s<0),1,0)-1 },
    `6` = { s1<-X[,"S1"];s2<-X[,"S2"]; 2*ifelse((s1>=0&s2>=-0.67)|(s1<0&s2<-0.67),1,0)-1 }
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

    fo <- train_rmst_cavboost(tr,tr$time,tr$status,tau,eta=0.05,max_depth=3,nr=50)
    po <- pred_subgroup(fo,te)

    ps <- tryCatch(predict(coxph(Surv(time,status)~.,data=tr[,c("time","status",f)]),type="lp"),error=function(e)NULL)
    if(!is.null(ps)) {
      st_ <- as.numeric(cut(ps,quantile(ps,seq(0,1,0.25),na.rm=T),include.lowest=TRUE))
      fs <- tryCatch(train_stratified_cavboost(tr,tr$time,tr$status,tau,stratum=st_,eta=0.1,max_depth=2,nr=50) ,error=function(e)NULL)
    } else fs <- NULL
    ps_te <- if(!is.null(fs)) pred_stratified(fs,te) else rep(0.5,nrow(te))

    ctrl <- tr[tr$trt01p==0,]; trt_d <- tr[tr$trt01p==1,]
    fm <- paste("Surv(time,status)~",paste(f,collapse="+"))
    fc <- suppressWarnings(tryCatch(coxph(as.formula(fm),data=ctrl,x=TRUE),error=function(e)NULL))
    ft <- suppressWarnings(tryCatch(coxph(as.formula(fm),data=trt_d,x=TRUE),error=function(e)NULL))
    vt <- if(!is.null(fc)&&!is.null(ft)) predict(ft,te,type="lp")-predict(fc,te,type="lp") else rep(0,nrow(te))

    cmf <- function(p,l) { r<-suppressMessages(roc(l,p,quiet=T))
      c(AUC=as.numeric(auc(r)),Acc=mean((p>.5)==l),
        FPR=if(sum(!l)>0)sum(p>.5&!l)/sum(!l)else NA,FNR=if(sum(l)>0)sum(p<=.5&l)/sum(l)else NA) }
    res <- rbind(res, data.frame(sc=label,rep=rep,method="OrigCB",t(cmf(po,l_te))))
    res <- rbind(res, data.frame(sc=label,rep=rep,method="StrCB",t(cmf(ps_te,l_te))))
    res <- rbind(res, data.frame(sc=label,rep=rep,method="VT",AUC=suppressMessages(roc(l_te,vt,quiet=T)$auc),Acc=NA,FPR=NA,FNR=NA))
  }, error=function(e) cat(sprintf("skip %s %d\n",label,rep)))
  cat(sprintf("%s done (%d)\n",label,sum(!is.na(res$AUC))))
  res
}

all <- data.frame()
for (sc in list(c(1,"S1_Linear"),c(3,"S3_U"),c(4,"S4_Enclave"),c(5,"S5_S"))) {
  all <- rbind(all, run_sc(as.numeric(sc[1]),sc[2]))
}
cat("\n=== Zhang scenarios (Setting 2, n_sim=30) ===\n")
for (sc in c("S1_Linear","S3_U","S4_Enclave","S5_S")) {
  cat(sprintf("\n--- %s ---\n",sc))
  cat(sprintf("%-12s %8s %8s %8s %8s\n","Method","AUC","Acc","FPR","FNR"))
  for(mn in c("OrigCB","StrCB","VT")){
    s<-all[all$sc==sc&all$method==mn,]
    cat(sprintf("%-12s %8.4f %8.4f %8.4f %8.4f\n",mn,
        mean(s$AUC,na.rm=T),mean(s$Acc,na.rm=T),mean(s$FPR,na.rm=T),mean(s$FNR,na.rm=T)))
  }
}
saveRDS(all,"zhang_results.rds")
