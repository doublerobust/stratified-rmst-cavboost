# MoE-K: Adaptive Stratification Results Summary

**Evaluation:** 989 scenarios × 5 methods (Oracle, Gate, CV, K=4, VT), across 8 scenario families and 4 sample sizes.

---

## 1. Overall AUC

| Method | Mean | Median | SD | p10 | p90 |
|--------|------|--------|----|-----|-----|
| Oracle | 0.7886 | 0.8056 | 0.128 | 0.601 | 0.954 |
| Gate | 0.7499 | 0.7528 | 0.139 | 0.563 | 0.937 |
| CV | 0.7395 | 0.7442 | 0.146 | 0.535 | 0.934 |
| 4-strata RMST | 0.7260 | 0.7317 | 0.149 | 0.525 | 0.932 |
| Virtual Twins | 0.7264 | 0.7270 | 0.126 | 0.564 | 0.895 |

**Gate-Oracle gap:** Mean −0.039, Median −0.003

---

## 2. AUC by Sample Size

| n_train | N | Oracle | Gate | CV | K4 | VT |
|---------|---|--------|------|-----|------|-----|
| 200 | 209 | 0.7096 | 0.6696 | 0.6548 | 0.6311 | 0.6647 |
| 300 | 206 | 0.7430 | 0.6979 | 0.6819 | 0.6740 | 0.6918 |
| 500 | 381 | 0.8136 | 0.7734 | 0.7665 | 0.7503 | 0.7466 |
| 1000 | 193 | 0.8734 | 0.8464 | 0.8399 | 0.8340 | 0.7904 |

**Gap from Oracle:**

| n_train | Gate | CV | K4 | VT |
|---------|------|----|-----|-----|
| 200 | −0.040 | −0.055 | −0.081 | −0.045 |
| 300 | −0.045 | −0.060 | −0.069 | −0.051 |
| 500 | −0.041 | −0.047 | −0.065 | −0.067 |
| 1000 | −0.027 | −0.037 | −0.043 | −0.083 |

The gate outperforms CV at every sample size. Gap shrinks with larger n (−0.040 → −0.027). VT collapses at n = 1000 (−0.083).

---

## 3. AUC by Scenario Family

| Family | N | Oracle | Gate | CV | K4 | VT |
|--------|---|--------|------|-----|------|------|
| additive | 128 | 0.764 | 0.734 | 0.719 | 0.724 | 0.707 |
| bump | 123 | 0.849 | 0.798 | 0.798 | 0.804 | 0.755 |
| cross | 134 | 0.747 | 0.680 | 0.680 | 0.688 | 0.643 |
| enclave | 137 | 0.793 | 0.767 | 0.755 | 0.694 | 0.765 |
| linear | 107 | 0.862 | 0.822 | 0.809 | 0.804 | 0.803 |
| radial | 134 | 0.740 | 0.719 | 0.706 | 0.654 | 0.718 |
| random | 99 | 0.793 | 0.751 | 0.743 | 0.722 | 0.748 |
| s_shaped | 125 | 0.780 | 0.747 | 0.726 | 0.735 | 0.692 |

**Gap from Oracle:**

| Family | Gate | CV | K4 | VT |
|--------|------|-----|------|-----|
| additive | −0.031 | −0.046 | −0.041 | −0.058 |
| bump | −0.051 | −0.052 | −0.046 | −0.095 |
| cross | **−0.067** | **−0.068** | −0.060 | **−0.103** |
| enclave | −0.026 | −0.039 | −0.094 | −0.028 |
| linear | −0.040 | −0.053 | −0.060 | −0.059 |
| radial | **−0.021** | −0.034 | −0.090 | **−0.022** |
| random | −0.042 | −0.051 | −0.076 | −0.045 |
| s_shaped | −0.035 | −0.057 | −0.055 | −0.088 |

**Best families for gate:** radial (−0.021), enclave (−0.026), additive (−0.031)
**Worst families:** cross (−0.067), bump (−0.051), random (−0.042)

---

## 4. AUC by Family & Sample Size

### Gate-Oracle gap

| Family | n=200 | n=300 | n=500 | n=1000 |
|--------|-------|-------|-------|--------|
| additive | −0.041 | −0.040 | −0.023 | −0.026 |
| bump | −0.032 | −0.086 | −0.059 | −0.017 |
| cross | −0.067 | −0.070 | **−0.070** | −0.058 |
| enclave | −0.025 | −0.043 | −0.026 | −0.013 |
| linear | −0.044 | −0.033 | −0.048 | −0.027 |
| radial | −0.019 | −0.016 | −0.019 | −0.034 |
| random | −0.060 | −0.039 | −0.033 | −0.039 |
| s_shaped | −0.031 | −0.045 | −0.050 | −0.008 |

### CV-Oracle gap

| Family | n=200 | n=300 | n=500 | n=1000 |
|--------|-------|-------|-------|--------|
| additive | −0.057 | −0.049 | −0.047 | −0.031 |
| bump | −0.065 | −0.076 | −0.040 | −0.037 |
| cross | −0.065 | −0.058 | −0.068 | −0.082 |
| enclave | −0.059 | −0.040 | −0.036 | −0.023 |
| linear | −0.053 | −0.057 | −0.062 | −0.024 |
| radial | −0.030 | −0.051 | −0.029 | −0.026 |
| random | −0.044 | −0.078 | −0.048 | −0.032 |
| s_shaped | −0.061 | −0.076 | −0.047 | −0.045 |

---

## 5. Gate vs CV: Head-to-Head

| Outcome | Count | % |
|---------|-------|---|
| Gate better (Δ > 0.005) | 285 | 28.8% |
| CV better (Δ < −0.005) | 197 | 19.9% |
| Tie (Δ within 0.005) | 507 | 51.4% |

Gate outperforms CV by ≈9 percentage points overall. In 129 more scenarios, CV wins marginally (within 0.005).

---

## 6. Large Failures

15% of scenarios (148/989) have gate lagging oracle by >0.10 AUC. Worst cases:

| Family | n | Oracle AUC | Gate AUC | Gap |
|--------|---|------------|----------|------|
| cross | 500 | 0.866 | 0.491 | −0.375 |
| s_shaped | 500 | 0.962 | 0.614 | −0.348 |
| bump | 300 | 0.929 | 0.612 | −0.317 |
| cross | 1000 | 0.932 | 0.639 | −0.294 |
| random | 500 | 0.428 | 0.139 | −0.290 |

Pattern: When oracle needed specific methods (m1–m4), gate fell back to method 6 (default CV), missing the optimal choice.

---

## 7. Feature Dictionary

All features are computed in one pass on training data — no cross-validation needed.

### A. Prognostic Signal
| Feature | Description |
|---------|-------------|
| `c_index` | Concordance index of null Cox model (covariates only, no treatment) |
| `score_skew` | Skewness of the linear predictor distribution |
| `score_kurt` | Kurtosis of the linear predictor distribution |
| `score_var` | Variance of the linear predictor |
| `score_q90_q10` | Ratio of 90th to 10th percentile of linear predictor |

### B. Interaction Structure
| Feature | Description |
|---------|-------------|
| `delta_c_index` | Improvement in C-index when adding treatment to null model |
| `prop_interact_sig` | Proportion of covariates with nominally significant (p<0.1) treatment interaction |
| `trt_main_p` | p-value of treatment main effect in Cox model |

### C. Treatment Effect Profile
| Feature | Description |
|---------|-------------|
| `te_bin_var` | Variance of RMST difference across prognostic bins (K=4) |
| `te_slope` | Slope of RMST difference vs. prognostic score (linear trend) |
| `te_quadratic_p` | p-value for quadratic deviation from linear TE trend |
| `te_max_diff` | Maximum pairwise difference in RMST across bins |
| `te_bin_event_rate_range` | Range of event rates across prognostic bins |

### D. Data Quality
| Feature | Description |
|---------|-------------|
| `censoring_rate` | Proportion of censored observations |
| `event_rate` | Proportion of events observed |
| `event_count_trt` | Event count in treatment arm |
| `event_count_ctrl` | Event count in control arm |
| `median_followup` | Median follow-up time |
| `mean_followup` | Mean follow-up time |
| `n` | Number of training samples |
| `p` | Number of covariates |
| `e_per_p` | Events per covariate (event_rate × n / p) |

### E. Sample Efficiency
| Feature | Description |
|---------|-------------|
| `bootstrap_ci_sd` | SD of bootstrap C-index estimates (50 replicates) |
| `bootstrap_ci_mean` | Mean of bootstrap C-index estimates |

### F. Model Behavior
| Feature | Description |
|---------|-------------|
| `orig_pred_var` | Variance of RMSTBoost treatment effect predictions |
| `orig_pred_mean` | Mean of RMSTBoost treatment effect predictions |
| `orig_ambiguity` | Proportion of predictions in [0.4, 0.6] (indeterminate zone) |

### H. Within-Arm Features
| Feature | Description |
|---------|-------------|
| `c_index_trt` | C-index within the treatment arm only |
| `c_index_ctrl` | C-index within the control arm only |
| `c_index_ratio` | Ratio of treatment to control C-index |

### I. TE Interaction Signals
| Feature | Description |
|---------|-------------|
| `te_int_max_z` | Maximum absolute z-score for treatment-covariate interactions |
| `te_int_mean_z` | Mean absolute z-score across all interactions |
| `te_int_prop_sig` | Proportion of interactions with \|z\| > 1.96 |

### J. Correlation Structure
| Feature | Description |
|---------|-------------|
| `corr_mean` | Mean absolute pairwise correlation among covariates |
| `corr_max` | Maximum absolute pairwise correlation |
| `corr_prop_high` | Proportion of pairwise correlations exceeding 0.5 |

---

## 8. Gate Model Architecture

- **Model:** ranger random forest (500 trees, impurity importance, probability output)
- **Training features:** 26 dataset-level meta-features (no oracle information)
- **Prediction targets:** Optimal K ∈ {1, 2, 3, 4, 5, VT} (CV-selected K)
- **Label source:** Oracle optimal K from simulation ground truth (training only)
- **Features explicitly excluded from training:** `auc_K1`–`auc_K5`, `optimal_K`, `oracle_auc`, `family`, `n_train`
