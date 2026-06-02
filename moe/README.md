# MoE-K: Adaptive Stratification via Meta-Learning Gate

See the [repo-level README](../README.md) for the full pipeline overview and results.

Key documents:
- [`draft-summary.md`](results/draft-summary.md) — Methods, results, statistical significance
- [`evaluate_gate.R`](evaluate_gate.R) — Ranger RF gate training + CV comparison
- [`gate_importance_viz.R`](gate_importance_viz.R) — Feature importance and activation heatmaps
- [`gate_features.R`](gate_features.R) — 39 meta-feature extractors
- [`results/gate_summary.pdf`](results/gate_summary.pdf) — AUC and match rate tables
- [`results/gate_vs_cv_histograms.png`](results/gate_vs_cv_histograms.png) — Distribution of AUC differences

### Quickstart (after simulation RDS files exist)

```bash
# 1. Extract training data from RDS
Rscript moe/extract_gate_data.R

# 2. Train gate + evaluate vs CV (10 parallel cores)
py moe/run_parallel_gate.py

# 3. Generate visualizations
Rscript moe/gate_importance_viz.R
```

### Running on Omen (Windows)

Requires Docker Desktop. See `run_on_omen_docker.sh` and `Dockerfile`.

```bash
docker build -t moe-k-sim -f moe/Dockerfile .
docker run --rm -v "$PWD/moe/results:/app/moe/results" --cpus 10 moe-k-sim \
    Rscript moe/evaluate_gate.R 0 10
```

### Model

The gate is a **ranger random forest** (500 trees, impurity importance, probability output) trained on 39 dataset-level features across 8 domains:

- Prognostic signal
- Treatment-effect heterogeneity
- Interaction structure
- Data maturity
- Sample regime
- Within-arm C-index
- Correlation structure
- Bootstrap uncertainty

The gate predicts the optimal K among {1, 2, 3, 4, 5, VT} directly from the training data, replacing 5-fold cross-validation.
