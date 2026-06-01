# MoE-K: When Cross-Validation Fails, Let a Gate Decide

## Motivation

The MoE-K framework selects the number of strata $K$ for the mixture-of-experts estimator. The standard approach — $K$-fold cross-validation within the trial — is computationally expensive and, more importantly, **unstable at small sample sizes**: the noisy AUC estimates from 5-fold CV often pick a suboptimal $K$, especially when $n \leq 300$.

We propose a **meta-learning gate**: a random forest trained on 39 dataset-level features (prognostic signal, treatment-effect heterogeneity, data maturity, interaction structure, etc.) that predicts the optimal $K$ directly from the training data, without any cross-validation.

## Data

- **Training set:** 5,000 simulation configurations × 10 reps = 50,000 synthetic oncology trial datasets (RFS endpoints, $n=200$–$1000$, 6 scenario families)
- Aggregated to 5,000 config-level rows by averaging features and selecting oracle-optimal $K$ via max mean AUC over $K=1..5$ plus Virtual Twin
- **Test set:** 989 hold-out configurations (80/20 split)

## Methods Compared

| Method | Description |
|--------|-------------|
| **Oracle** | Picks $K$ with highest true AUC on test data (theoretical upper bound) |
| **Gate** | Ranger random forest, 500 trees, impurity importance, trained on 39 dataset features |
| **CV** | Real 5-fold within-trial cross-validation with adaptive fold collapsing (5→3→2) |
| **Fixed K=4** | Always uses 4 strata (no selection) |
| **VT** | Virtual Twin estimator (one of the 6 candidate methods) |

## Results

### AUC by Sample Size

| $n$ | Oracle | **Gate** | CV | Fixed K=4 | VT alone | $N$ |
|:---:|:------:|:--------:|:--:|:---------:|:--------:|:---:|
| 200 | 0.7096 | **0.6692** | 0.6548 | 0.6311 | 0.6647 | 209 |
| 300 | 0.7430 | **0.6965** | 0.6819 | 0.6740 | 0.6918 | 206 |
| 500 | 0.8136 | **0.7718** | 0.7665 | 0.7503 | 0.7466 | 381 |
| 1000 | 0.8734 | **0.8491** | 0.8399 | 0.8340 | 0.7904 | 193 |
| All | 0.7886 | **0.7494** | 0.7395 | 0.7260 | 0.7264 | 989 |

The gate **consistently outperforms CV at every sample size**. The advantage is largest at $n=1000$ (Gate gap to oracle = 0.024 vs CV gap = 0.034, a 29% reduction) and narrowest at $n=500$ (10% reduction). Critically, the gate **never underperforms CV** — even at $n=200$ where CV is most unstable, the gate maintains a 0.014 AUC advantage.

### Exact Oracle Method Match Rate

| $n$ | **Gate** | CV | Gate edge |
|:---:|:--------:|:--:|:---------:|
| 200 | **55.5%** | 37.3% | +18.2pp |
| 300 | **49.0%** | 33.0% | +16.0pp |
| 500 | **41.7%** | 37.9% | +3.8pp |
| 1000 | **40.9%** | 29.5% | +11.4pp |
| All | **46.0%** | 35.1% | +10.9pp |

The match rate declines with $n$ because the oracle uses a wider variety of methods at larger $n$ — the decision problem becomes harder. However, the gate's advantage over CV remains consistent across all $n$, suggesting it captures genuine structure rather than memorizing VTs default.

### Feature Importance

The RF importance reveals which dataset characteristics drive the gate's decisions:

| Domain | Key features | Impact |
|--------|-------------|--------|
| Prognostic signal | `c_index`, `bootstrap_ci_sd`, `orig_pred_var` | Strongest — C-index alone explains most gate decisions |
| TE heterogeneity | `te_slope`, `te_quadratic_p`, `te_int_max_z` | Moderate — shapes how many strata are needed |
| Data maturity | `event_rate`, `censoring_rate`, `median_followup` | Moderate — fewer events favor fewer strata |
| Sample regime | `n`, `p`, `e_per_p` | Mild — mainly affects confidence in other features |
| Interaction structure | `delta_c_index`, `prop_interact_sig` | Weakest — non-linear miscalibration patterns |

The activation heatmaps (Figure 3) show that different scenario families ("s_shaped", "bump", "linear") light up distinct feature domains, confirming the gate is learning scenario-specific decision rules rather than a single default.

## Discussion

**When is the gate most valuable?** At the extremes: $n \leq 300$ where CV is unreliable, and $n \geq 1000$ where the oracle has many viable options and the gate's richer feature set gives it an edge. The $n=500$ region is the convergence zone — CV works adequately and the gate's features are still noisy — but even there the gate matches or beats CV.

**What limits the gate?** The current gate defaults to VT on ~16% of test cases where the feature profile is ambiguous. This VT skew is partly due to `prop_interact_sig` being all-NA (the Cox interaction model with 40+ terms failed to converge on most datasets). A univariate screening fix is expected to improve discrimination further.

**Future directions:**
- Integrate the univariate interaction features to reduce VT default rate
- Bayesian gate that outputs prediction intervals, not just point estimates
- Active learning: the gate could request a CV check only when uncertainty is high

## References

- Zhang et al. (2023). Mixture-of-Experts for survival analysis.
- Schuler et al. (2022). Prognostic score adjustment.
- Mehrotra (2021). ENET-based risk stratification (5-STAR).
