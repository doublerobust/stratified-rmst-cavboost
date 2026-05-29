# Stratified RMST CAVBoost

Stratified RMST value function boosting for subgroup identification in
time-to-event outcomes.  Augments the original CAVBoost gradient by
computing the value function within prognostic strata, blocking
prognostic confounding and improving subgroup identification accuracy.

## Key Formula

Value function:

$$V_{\text{strat}} = \sum_{k=1}^K W_k\, d_k^{(1)} \;-\; \sum_{k=1}^K (n_k - W_k)\, d_k^{(2)}$$

Gradient (XGBoost minimizes):

$$g_j = -p_j(1-p_j)\Big[d_k^{(1)} + d_k^{(2)} + W_k\frac{\partial d_k^{(1)}}{\partial p_j} - (n_k-W_k)\frac{\partial d_k^{(2)}}{\partial p_j}\Big]$$

## Repository Structure

```
├── R/                          # Implementation
│   ├── stratified_cavboost.R   # Main gradient implementation
│   └── test_stratified_gradient.R  # Gradient verification
├── sim/                        # Simulation scripts
│   ├── run_holdout.R           # Primary DGP comparison
│   ├── run_vt.R                # Virtual twin comparison
│   ├── run_sb_prog.R           # SubgroupBoost-style DGP
│   ├── run_zhang.R             # Zhang complex boundary scenarios
│   └── compare_prog_scores.R   # Prognostic score method comparison
├── methodology/                # Writeup
│   ├── stratified-rmst-boosting.pdf
│   └── stratified-rmst-boosting.tex
└── README.md
```

## Key Results (hold-out AUC)

| Scenario | Original | Stratified | VT | Gain |
|----------|:--------:|:----------:|:--:|:----:|
| *Primary DGP* | 0.602 | **0.635** | 0.705 | +3.3 pts |
| *Pure predictive* | 0.620 | **0.656** | 0.706 | +3.6 pts |
| *SubgroupBoost DGP* | 0.679 | **0.740** | 0.762 | +6.1 pts |
| S1: Linear | 0.798 | **0.964** | 0.956 | +16.6 pts |
| S3: U-shaped | 0.918 | **0.977** | 0.509 | +5.9 pts |
| S4: Enclave | 0.844 | 0.765 | 0.506 | -7.9 pts |
| S5: S-shaped | 0.918 | 0.849 | 0.596 | -6.9 pts |

## Dependencies

- R 4.6+ with: survival, xgboost, mvtnorm, pROC, pseudo
- For the original CAVBoost reference: `rmst_cavboost_clean.R` from
  the [CAVBoost](https://github.com/sambiostat/CAVBoost) repository.

## Reference

Zhang, P., Liu, P., Chen, X., Ma, J., & Shentu, Y. (2020).
*A nonparametric method for value function guided subgroup identification
via gradient tree boosting for censored survival data.*
Statistics in Medicine, 39(28), 4133--4146. PMID: 32786155.
