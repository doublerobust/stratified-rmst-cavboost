# Dependencies & Onboarding

## System Requirements
- R >= 4.2
- 16 GB+ RAM recommended (500 obs × 52 covariates × 200-tree ranger)

## Required R Packages
```r
install.packages(c(
  "xgboost",       # gradient tree boosting
  "glmnet",        # elastic net Cox for prognostic strata
  "ranger",        # survival random forests for virtual twin
  "survival",      # Surv objects, Cox PH
  "mvtnorm",       # multivariate normal for simulation DGP
  "pROC",          # AUC computation (optional — manual AUC fallback exists)
  "dplyr"          # data manipulation
))
```

## File Structure
```
stratified-rmst-cavboost/
├── R/
│   ├── stratified_cavboost.R    # main implementation
│   └── rmst_cavboost_clean.R    # original CAVBoost (baseline)
├── run_one_rep_v2.R             # single-rep runner (fresh R process each call)
├── run_main_comparison_50.R     # in-process 50-rep runner
├── results_final/               # per-rep RDS outputs (generated)
└── specs/
    └── foundation-moe-subgroup.md
```

## Running the Simulation

### Single rep (testing):
```bash
Rscript run_one_rep_v2.R <scenario 1-6> <rep 1-50> <output_path>
```

### Full 50-rep × 6-scenario (sequential, about 3h):
```bash
for sc in 1 2 3 4 5 6; do
  for rep in $(seq 1 50); do
    Rscript run_one_rep_v2.R $sc $rep results_final/sc${sc}_rep${rep}.rds
  done
done
```

### Parallel (on Omen RTX 5090):
```bash
# Use GNU parallel to run 8+ reps simultaneously
parallel -j 8 Rscript run_one_rep_v2.R {1} {2} results_final/sc{1}_rep{2}.rds ::: \
  1 2 3 4 5 6 ::: $(seq 1 50)
```

## Key Parameter Agreements
- **Prognostic strata**: 5-fold cross-fitted `cv.glmnet(family="cox", alpha=0.5)`, all 52 covariates
- **Original CAVBoost**: `eta=0.05, max_depth=3, nrounds=50`
- **Stratified CAVBoost**: `eta=0.1, max_depth=2, nrounds=50`
- **Virtual Twin**: `ranger(num.trees=200, min.node.size=10)`, per-arm survival forests, RMST via trapezoidal integration on τ=30
- **Zhang scenarios**: `n_train=500, n_test=2000, τ=30, rho=1/3`, Setting 2 prognostic strength (β=0.4)

## Reproducing the Paper Results
```r
# Collect results from per-rep RDS files
files <- list.files("results_final", "\\.rds$", full.names = TRUE)
res <- do.call(rbind, lapply(files, readRDS))
aggregate(cbind(orig_auc, strat_auc, vt_auc) ~ sc, data = res, mean)
```
