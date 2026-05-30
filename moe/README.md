# MoE-K: Adaptive Stratification for RMST Boosting

## The Core Idea

Stratified RMST boosting has one tuning parameter: **K**, the number of
prognostic strata.  K controls the bias-variance tradeoff:

- **K = 1**: Original RMSTBoost — no stratification. High bias (prognostic
  confounding can dilute the gradient), low variance (all data pooled).
- **K = 4** (default): Moderate stratification — typical choice for
  n = 500. Reduces prognostic confounding at the cost of smaller per-stratum
  samples.
- **K = 6 or 8**: Strong stratification — lower bias but each stratum may
  have too few events for stable KM estimation, especially under heavy
  censoring.

The optimal K depends on dataset properties that are **observable at
training time** — prognostic signal strength, interaction complexity,
sample size, censoring rate, and the shape of the treatment-effect
boundary.  Rather than fixing K = 4, we learn a **gating function** that
predicts the optimal K from these observable features.

## Advantages Over 3-Expert MoE

| | 3-Expert MoE | Within-Family K-Selection |
|---|---|---|
| Model class | Orig, StratCF, VT (heterogeneous) | RMSTBoost only (homogeneous) |
| Gate's job | "Which model class?" | "How much stratification?" |
| Failure mode | Wrong expert → disaster | Suboptimal K → slight efficiency loss |
| Regulatory path | Hard (multiple models) | Easier (one tunable method) |
| Computational cost | 3x model fitting | 1x model fit at chosen K |

## Feature Library (~40 dataset-level descriptors)

All computable in one pass on the training data (no CV):

### A. Prognostic Signal (fit Cox on Z only, no A)
1. C-index of Cox model
2. Skewness of linear predictor (prognostic score)
3. Kurtosis of linear predictor
4. Dip test p-value (multimodality of score)
5. Variance of score across patients
6. Ratio of top 10% to bottom 10% hazard

### B. Interaction Structure (compare Z-only vs Z+A models)
7. ΔC-index: Cox(Z × A) minus Cox(Z only)
8. Log-rank test p-value: Cox residuals split by treatment
9. Interaction F-test p-value from Cox with Z + Z:A
10. Proportion of β_trt interactions with p < 0.1
11. Concordance of stratified KM (by prognostic quartile) curves

### C. Treatment Effect Profile (bin prognostic score)
12. Variance of RMST difference across 4 bins
13. Slope of RMST difference vs. mean-bin score
14. Deviation from linearity (quadratic contrast) in RMST bins
15. Maximum pairwise RMST diff between bins
16. Within-bin event rate range

### D. Data Quality
17. Censoring rate (overall and per arm)
18. Event count and event rate per arm
19. Mean / median follow-up time
20. Proportion of events past τ/2
21. Ratio of events to covariates (E/p)

### E. Sample Efficiency Indicators
22. Leave-one-out stability of Cox coefficients
23. Bootstrap variance of C-index estimate
24. Effective sample size (information-based)

### F. Oracle-Derived (from model outputs on validation split)
25. Orig RMSTBoost internal prediction variance
26. StratCF (K=4) internal prediction variance
27. Pairwise prediction MSE between Orig and StratCF
28. Mean subgroup probability (p_i) for each model
29. Proportion of patients with p_i in [0.4, 0.6] (ambiguity zone)

## Simulation Data Generator

Sample across continuous axes using a parametric boundary generator:

| Parameter | Range | Notes |
|---|---|---|
| Boundary family | {linear, additive, bump, enclave, S-shaped, cross, radial} | Plus random combinations |
| # predictive variables | 1–10 | True treatment-effect modifiers |
| # prognostic variables | 0–52 | Affect baseline hazard only |
| Overlap of predictive/prognostic | {none, partial, complete} | |
| Prognostic strength (b₀) | 0–2 | Baseline log-HR |
| Prognostic form | {linear, quadratic, step, saturating} | |
| Sample size (train) | 100, 200, 300, 500, 1000, 2000 | |
| Censoring rate | 10%–70% | At τ = 30 |
| Covariate correlation | {low, moderate, high} | Block-diagonal or AR(1) |

**Target:** 200 diverse configurations × 5 reps = 1,000 runs.

## Gate Training

### Architecture
Regularized multinomial logistic regression (LASSO) with softmax output
over K ∈ {1, 2, 3, 4, 5}:

```
P(K = k | x) = exp(β_k · x) / Σ_j exp(β_j · x)
```

where x is the 40-dimensional feature vector.

### Training Labels
For each simulation rep:
1. Fit RMSTBoost for each K ∈ {1, 2, 3, 4, 5} on training data
2. Evaluate AUC on held-out test data
3. Label = argmax_k AUC_k

### Optimization
- Loss: categorical cross-entropy
- Regularization: 10-fold CV to select LASSO λ
- Class balancing: K=1 and K=5 may be rare; use synthetic oversampling or
  weighted loss

### Validation
1. **Within-family**: 80/20 train/test split across all reps
2. **Cross-boundary**: train on linear+additive+bump families, test on
   enclave+S-shaped+cross (hard generalization)
3. **Null simulation**: train gate on active-treatment data; test on data
   with no treatment effect → gate should output K ≈ 1 (no stratification
   needed when no subgroups exist)

## Expected Deliverables
1. Parametric scenario sampler (R)
2. Gate feature extraction pipeline (R)
3. Trained LASSO gate (small coefficient vector, ~5–15 non-zero)
4. Validation report showing:
   - AUC gain from adaptive K vs. fixed K=4
   - Which features are selected by LASSO (interpretability)
   - Gate performance on held-out boundary families

## File Structure (on moe-integration branch)
```
moe/
├── scenario_generator.R    # Parametric scenario sampler
├── moe_simulation.R        # Full simulation (data gen → fit → eval → feature extract)
├── gate_features.R          # Feature computation from training data
├── train_gate.R            # LASSO gate training
├── evaluate_gate.R         # Validation on held-out scenarios
├── gate_coefficients.csv   # Trained gate output
└── README.md               # This plan
```

## Open Questions
1. How much does adaptive K gain over fixed K=4? (Primary result.)
2. Which 5–10 features does LASSO select? (Interpretability finding.)
3. Does the gate generalize to boundary families it never saw? (Robustness.)
4. can we run 1,000 simulation runs in ~1 hour on Omen? (Feasibility.)
5. Should K be an integer (discrete choice) or a real-valued smoothing
   parameter? (Design decision.)
