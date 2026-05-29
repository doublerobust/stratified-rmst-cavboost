# Run ONE rep of holdout selection, save result, exit.
# Usage: Rscript run_one_rep.R <scenario> <rep> <outfile>

args <- commandArgs(trailingOnly = TRUE)
sc <- as.numeric(args[1])
rep <- as.numeric(args[2])
outfile <- args[3]

library(pROC); library(mvtnorm)
source("R/stratified_cavboost.R")
source("~/.openclaw/workspace/CAVBoost/rmst_cavboost_clean.R")

tau <- 30; n_tr <- 500; n_te <- 2000; nr <- 50
Smat <- matrix(1/3, 52, 52); diag(Smat) <- 1
feat_prog <- c("Z1","Z2","Z3","Z4","S1","S2")

set.seed(42 + rep * 1000 + sc * 100)
X <- rmvnorm(n_tr + n_te, sigma = Smat)
colnames(X) <- c(paste0("z",1:50),"S1","S2")
colnames(X)[1:4] <- c("Z1","Z2","Z3","Z4")
A <- rbinom(n_tr + n_te, 1, 0.5); b0 <- sqrt(6); s0 <- 0.4
Zb <- X[,1:4,drop=F] %*% c(0.4,0.4,0.4,0.4)

if (sc == 4) {
  s1 <- X[,"S1"]
  te <- 2*((-1.07 <= s1 & s1 < 1.07) & (-1.07 <= X[,"S2"] & X[,"S2"] < 1.07)) - 1
} else {
  s <- X[,"S1"]
  te <- 2*ifelse(s >= 0.67 | (-0.67 <= s & s < 0), 1, 0) - 1
}

Zb <- -Zb^2
T <- exp(b0 + A * te + Zb + s0 * rnorm(n_tr + n_te))
C <- pmin(30, rexp(n_tr + n_te, rate = -log(0.9)/12))
U <- pmin(T, C, tau); st <- as.numeric(T <= C)
oracle <- as.numeric(te > 0)

if (length(unique(oracle)) < 2) { saveRDS(NULL, outfile); q() }

tr <- data.frame(X[1:n_tr,], trt01p = A[1:n_tr], time = U[1:n_tr], status = st[1:n_tr])
te_df <- data.frame(X[-(1:n_tr),], trt01p = A[-(1:n_tr)], time = U[-(1:n_tr)], status = st[-(1:n_tr)])
l <- oracle[-(1:n_tr)]

auc_ <- function(p,l) tryCatch(as.numeric(pROC::auc(pROC::roc(l,p,quiet=TRUE,direction="<"))), error=function(e) NA)

# Original
fo <- tryCatch(train_rmst_cavboost(tr,tr$time,tr$status,tau,eta=0.05,max_depth=3,nr=nr), error=function(e) NULL)
po <- if(!is.null(fo)) pred_subgroup(fo, te_df) else rep(0.5, nrow(te_df))
rm(fo)

# Stratified (cross-fitted)
st_ <- tryCatch(crossfit_prognostic_strata(tr, feat_prog, nfold=5, K=4, seed=rep*100+sc), error=function(e) NULL)
fs <- if(!is.null(st_)) tryCatch(train_stratified_cavboost(tr,tr$time,tr$status,tau,stratum=st_,eta=0.1,max_depth=2,nr=nr), error=function(e) NULL) else NULL
ps <- if(!is.null(fs)) pred_stratified(fs, te_df) else rep(0.5, nrow(te_df))
rm(fs)

# Holdout selection
sel <- tryCatch(select_model_by_holdout(tr,"time","status",tau,stratum=st_,features=feat_prog), error=function(e) NULL)
psel <- if(!is.null(sel)) {
  if(sel$winner == "original") pred_subgroup(sel$fit, te_df) else pred_stratified(sel$fit, te_df)
} else rep(0.5, nrow(te_df))

winner <- if(!is.null(sel)) sel$winner else "none"
gain_o <- if(!is.null(sel)) sel$gain_orig else NA
gain_s <- if(!is.null(sel)) sel$gain_strat else NA

result <- data.frame(
  sc = sc, rep = rep,
  orig_auc = auc_(po, l), strat_auc = auc_(ps, l), select_auc = auc_(psel, l),
  winner = winner, gain_o = gain_o, gain_s = gain_s
)
saveRDS(result, outfile)
cat(sprintf("rep %d done (winner: %s)\n", rep, winner))
