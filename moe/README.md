# MoE-K: Adaptive Stratification via Meta-Learning Gate

Pre-trained meta-learning gate that selects the optimal number of stratification strata K for RMST boosting, replacing per-dataset 5-fold cross-validation. The gate operates on 26 dataset-level summary statistics and is trained on 1,533 synthetic scenarios.

## Results Summary

Evaluation on 989 independent scenarios across 8 families and 4 sample sizes:

| Method | Mean AUC | Median AUC | Gap vs Oracle |
|--------|----------|------------|---------------|
| Oracle | 0.789 | 0.806 | — |
| **Gate** | **0.750** | **0.753** | **−0.039** |
| CV | 0.740 | 0.744 | −0.048 |
| 4-strata RMST | 0.726 | 0.732 | −0.063 |
| Virtual Twins | 0.726 | 0.727 | −0.063 |

**Gate outperforms CV at every sample size** (285 vs 197 head-to-head wins, 51.4% ties).

### By Sample Size

| n | Gate | CV | Gate-Over-CV |
|---|------|----|--------------|
| 200 | 0.670 | 0.655 | **+0.015** |
| 300 | 0.698 | 0.682 | **+0.016** |
| 500 | 0.773 | 0.767 | **+0.007** |
| 1000 | 0.846 | 0.840 | **+0.007** |

### Key Findings

- **Gate beats CV** by ≈9 percentage points in head-to-head comparisons (29% vs 20%)
- **Gap to oracle shrinks** with larger n: from −0.040 at n=200 to −0.027 at n=1000
- **Best families:** radial (−0.021), enclave (−0.026), additive (−0.031)
- **Worst families:** cross (−0.067), bump (−0.051), random (−0.042)
- **15% of scenarios** have gate lagging oracle by >0.10 AUC (usually when optimal method is m1–m4)

## Documents

| File | Description |
|------|-------------|
| [`results/draft-summary.md`](results/draft-summary.md) | Full results with tables by family, n, and head-to-head comparisons |
| [`results/gate_vs_cv_violin.png`](results/gate_vs_cv_violin.png) | Violin plots comparing Gate vs CV AUC distributions |
| [`../introduction_draft.md`](../introduction_draft.md) | Manuscript introduction and framing |
| [`evaluate_gate.R`](evaluate_gate.R) | Ranger RF gate training + CV comparison |
| [`gate_importance_viz.R`](gate_importance_viz.R) | Feature importance and activation heatmaps |
| [`gate_features.R`](gate_features.R) | 26 meta-feature extractors across 9 domains |

## Quickstart

```bash
# 1. Extract training data from RDS
Rscript moe/extract_gate_data.R

# 2. Train gate + evaluate vs CV (10 parallel cores)
py moe/run_parallel_gate.py

# 3. Generate visualizations
Rscript moe/gate_importance_viz.R
Rscript moe/gate_violin_plots.R
```

## Model

The gate is a **ranger random forest** (500 trees, probability output) trained on 26 dataset-level features across 9 domains:

- Prognostic signal (c_index, score skewness/kurtosis/variance)
- Interaction structure (ΔC-index, proportion significant interactions)
- Treatment-effect profile (TE variance across bins, TE slope, quadratic p)
- Data quality (censoring rate, event counts, e/p ratio)
- Bootstrap uncertainty (C-index bootstrap SD)
- Within-arm features (treatment/control C-index, ratio)
- TE interaction signals (max/mean z-scores, proportion significant)
- Correlation structure (mean/max correlation, proportion high)
- Original model behavior (prediction variance, mean, ambiguity)

**No oracle information leaks into features.** All 26 features are computed in one pass on the training data without knowledge of the true treatment effect. The `auc_K1`–`auc_K5` and `optimal_K` columns are used only as evaluation targets, not as predictors.

## Feature Explanations

Detailed descriptions of all 26 features are in the [results/draft-summary.md §7](results/draft-summary.md#7-feature-dictionary).

## Running on Omen (Windows)

Requires Docker Desktop. See `run_on_omen_docker.sh` and `Dockerfile`.

```bash
docker build -t moe-k-sim -f moe/Dockerfile .
docker run --rm -v "$PWD/moe/results:/app/moe/results" --cpus 10 moe-k-sim \
    Rscript moe/evaluate_gate.R 0 10
```
