# Stratified RMST CAVBoost

Stratified RMST value function boosting for subgroup identification in
time-to-event outcomes.  Augments the original CAVBoost gradient by
computing the value function within prognostic strata, blocking
prognostic confounding and improving subgroup identification accuracy.

The prognostic score is estimated via **5-fold cross-fitted elastic net
Cox** (`cv.glmnet`, `alpha=0.5`) using **all 52 covariates**.  No
collapsing or variable subsetting — experience shows that using fewer
covariates silently changes the score and produces misleading comparisons.

## Key Formula

Value function:

$$V_{\text{strat}} = \sum_{k=1}^K W_k\, d_k^{(1)} \;-\; \sum_{k=1}^K (n_k - W_k)\, d_k^{(2)}$$

Gradient (XGBoost minimizes):

$$g_j = -p_j(1-p_j)\Big[d_k^{(1)} + d_k^{(2)} + W_k\frac{\partial d_k^{(1)}}{\partial p_j} - (n_k-W_k)\frac{\partial d_k^{(2)}}{\partial p_j}\Big]$$

## Repository Structure

```
├── R/
│   ├── stratified_cavboost.R      # Main gradient + crossfit implementation
│   ├── rmst_cavboost_clean.R      # Original CAVBoost (baseline comparator)
│   └── test_stratified_gradient.R # Gradient verification
├── methodology/
│   ├── stratified-rmst-boosting.tex  # Manuscript source
│   └── stratified-rmst-boosting.pdf  # Rendered manuscript
├── specs/
│   └── foundation-moe-subgroup.md    # Future work: Foundation MoE spec
├── run_main_comparison_50.R       # In-process 50-rep runner (Orig vs Strat vs VT)
├── run_clean_50.R                 # Clean sequential runner (CSV output, gc() between reps)
├── run_one_rep_v2.R               # Single-rep runner (fresh R process)
├── run_sequential.R               # Subprocess sequential runner (calls v2)
├── DEPENDENCIES.md                # Onboarding & installation guide
└── README.md
```

## Running the Simulation

### Single rep (testing):
```bash
Rscript run_one_rep_v2.R <scenario 1-6> <rep 1-50> <output_path>
```

### Full 50-rep × 6-scenario (sequential, about 3h):
```bash
Rscript run_clean_50.R > clean.log 2>&1 &
```

### Full run via subprocess:
```bash
Rscript run_sequential.R
```

### Parallel (on a multi-core machine):
```bash
parallel -j 8 Rscript run_one_rep_v2.R {1} {2} results/sc{1}_rep{2}.rds ::: \
  1 2 3 4 5 6 ::: $(seq 1 50)
```

### In-process (memory-intensive, faster):
```bash
Rscript run_main_comparison_50.R
```

## Key Parameter Agreements

| Parameter | Value |
|-----------|-------|
| Prognostic strata | 5-fold cross-fitted `cv.glmnet(family="cox", alpha=0.5)`, all 52 covariates |
| Original CAVBoost | `eta=0.05, max_depth=3, nrounds=50` |
| Stratified CAVBoost | `eta=0.1, max_depth=2, nrounds=50` |
| Virtual Twin | `ranger(num.trees=200, min.node.size=10)`, per-arm RF, RMST via trapezoidal integration on τ=30 |
| Zhang scenarios | `n_train=500, n_test=2000, τ=30, ρ=1/3`, Setting 2 prognostic strength |

## Dependencies

- R 4.2+ with: xgboost, glmnet, ranger, survival, mvtnorm, pROC
- See `DEPENDENCIES.md` for full installation instructions.

## Reference

Zhang, P., Liu, P., Chen, X., Ma, J., & Shentu, Y. (2020).
*A nonparametric method for value function guided subgroup identification
via gradient tree boosting for censored survival data.*
Statistics in Medicine, 39(28), 4133--4146. PMID: 32786155.
