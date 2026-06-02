# Stratified RMST CAVBoost

Stratified RMST value function boosting for subgroup identification in time-to-event outcomes — with a **Mixture-of-Experts (MoE-K) meta-learning gate** that automatically selects the optimal number of strata.

## MoE-K Gate: When Cross-Validation Fails

The standard approach to selecting the number of strata $K$ is within-trial cross-validation. But CV is noisy at small sample sizes and expensive at large ones. The **MoE-K gate** replaces CV with a random forest trained on 39 dataset-level meta-features that predicts the optimal $K$ directly.

### Key Results (989 test configurations)

| Method | Mean AUC | Exact match vs oracle |
|--------|:--------:|:--------------------:|
| Oracle (theoretical upper bound) | 0.7886 | — |
| **Gate (ranger RF)** | **0.7494** | **46.0%** |
| 5-fold CV | 0.7395 | 35.1% |
| Fixed K=4 | 0.7260 | — |
| VT alone | 0.7264 | — |

The gate beats CV at every sample size ($n=200$–$1000$), with a mean AUC advantage of +0.01 (paired t-test $p = 6.4 \times 10^{-5}$). Full results in `moe/results/draft-summary.md`.

### Feature Domains

The gate uses 39 features across 8 domains: prognostic signal, treatment-effect heterogeneity, interaction structure, data maturity, sample regime, within-arm C-index, correlation structure, and bootstrap uncertainty. The [global importance plot](https://github.com/doublerobust/stratified-rmst-cavboost/blob/moe-integration/moe/results/gate_global_importance.pdf) shows C-index and bootstrap CI width are the strongest drivers.

## Repository Structure

```
├── R/                           # CAVBoost estimators
│   ├── stratified_cavboost.R    # Main gradient + crossfit
│   ├── rmst_cavboost_clean.R    # Original CAVBoost (comparator)
│   └── test_stratified_gradient.R
├── moe/                         # MoE-K gate pipeline
│   ├── gate_features.R           # 39 meta-feature extractors
│   ├── evaluate_gate.R           # Gate training + CV comparison
│   ├── gate_importance_viz.R     # Feature importance + activation heatmaps
│   ├── scenario_generator.R      # Parametric dataset sampler
│   ├── extract_gate_data.R       # Aggregate RDS → training CSV
│   ├── run_parallel_gate.py      # 10-core parallel launcher
│   ├── run_parallel_simulation.py
│   ├── run_on_omen.sh            # Omen launcher (legacy)
│   ├── run_on_omen_docker.sh     # Omen launcher (Docker)
│   ├── Dockerfile                # Reproducible R 4.6.0 environment
│   ├── build.sh                  # Docker build + push
│   └── results/
│       ├── gate_training_data.csv        # Training data (5000 configs)
│       ├── gate_evaluation.csv           # Test results (989 splits)
│       ├── gate_summary.pdf              # AUC + match rate tables
│       ├── gate_vs_cv_histograms.png     # Distribution of AUC differences
│       ├── gate_global_importance.pdf    # Feature importance by domain
│       ├── gate_activation_heatmap.pdf   # Per-dataset activation map
│       ├── gate_brain_slice.pdf          # Single-dataset fMRI-style view
│       └── draft-summary.md              # Methods + results narrative
├── code/                        # Older analysis scripts
└── README.md
```

## Running the Gate Pipeline

### 1. Data extraction (from simulation RDS files)
```bash
Rscript moe/extract_gate_data.R
```

### 2. Train gate + run CV comparison
```bash
# Sequential:
Rscript moe/evaluate_gate.R

# Parallel (10 cores):
python3 moe/run_parallel_gate.py
```

### 3. Generate visualizations
```bash
Rscript moe/gate_importance_viz.R
```

### Running on Omen (Windows) — Docker

```bash
# One-time build
docker build -t moe-k-sim -f moe/Dockerfile .

# Run evaluation
docker run --rm -v "$PWD/moe/results:/app/moe/results" --cpus 10 moe-k-sim \
    Rscript moe/evaluate_gate.R 0 10
```

See `moe/run_on_omen_docker.sh` for the full launcher.

## Dependencies

- R 4.2+ with: ranger, glmnet, survival, xgboost, mvtnorm, data.table, ggplot2, gridExtra, RColorBrewer, viridisLite
- Python 3 with: pandas, numpy (for parallel launchers)
- Docker (optional, for Windows/Omen runs)

## Key Algorithm Parameters

| Parameter | Value |
|-----------|-------|
| Gate model | Ranger RF, 500 trees, impurity importance |
| Gate features | 39 meta-features from 8 domains |
| CV comparator | Real 5-fold with adaptive collapse (5→3→2) |
| Candidate K | 1, 2, 3, 4, 5, VT (Virtual Twin) |
| Base estimator | RMST CAVBoost (XGBoost, eta=0.05, max_depth=3) |
| Simulation reps | 5000 configs × 10 reps = 50,000 datasets |
| Sample sizes | n = 200, 300, 500, 1000 |
| Scenario families | s_shaped, bump, linear, radial, enclave, cross, additive, random |

## MoE-K Reference

For the underlying stratified RMST estimator and original RMST Boosting:

Zhang P, Ma J, Chen X, Shentu Y. A nonparametric method for value function guided subgroup identification via gradient tree boosting for censored survival data. Statistics in Medicine. 2020;39:4133–4146. https://doi.org/10.1002/sim.8714
