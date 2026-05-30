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

## MoE-K vs. Cross-Validation

The most natural baseline: fit RMSTBoost for K = 1..5, pick the one with
the best cross-validated AUC.  MoE-K should outperform CV because:

| | MoE-K | CV |
|---|---|---|
| **Training labels** | Oracle (true test AUC on 2,000 held-out patients) | Noisy CV estimates with large SE |
| **At n = 200** | Gate learns small-sample patterns | CV can't distinguish K=2 from K=4 (SE ≈ 0.02 on a 0.01 AUC gap) |
| **At n = 2000** | May be comparable | CV works well |
| **Cost (new dataset)** | One Cox fit + 40 summary stats + matrix multiply (~0.5s) | 25 RMSTBoost fits (~5 min) |
| **Interpretability** | ~5–15 non-zero LASSO coefficients | "CV said K=4" |

**Hypothesis:** MoE-K's advantage over CV is largest at small n and
diminishes as sample size grows.  The gate learns from oracle labels
across 1,000 simulation reps, effectively pooling information across
datasets in a way that a single-dataset CV cannot.

This is itself a publishable secondary result: **the sample size regime
where pre-trained gates beat cross-validation.**

## Training Data Design

### Sample Size

Train the gate on the same sample sizes we deploy on.  Early-phase
oncology trials typically run n = 200–1000.  Simulate at multiple n levels
and include **n as an explicit gate feature**, so the gate learns:

- "At n = 200, even strong prognostic signal is noisy → conservative K = 2–3"
- "At n = 1000, clean signal → aggressive K = 5"
- "At n = 500, moderate signal → K = 4 unless ΔC-index is very small"

Sampling plan (n × 5 reps = 1,000 runs total):

| n | Configurations | Reps | Total runs | Rationale |
|---|---|---|---|---|
| 200 | 40 | 5 | 200 | Smallest realistic, worst-case for CV |
| 300 | 40 | 5 | 200 | Typical Phase Ib/IIa |
| 500 | 80 | 5 | 400 | Our current simulation N, sweet spot |
| 1000 | 40 | 5 | 200 | Larger, CV starts working well |

The gate's LASSO will decide whether n is useful as a main effect,
whether it interacts with other features, or whether it's dominated by
other signals.

### Parametric Scenario Sampler

Sample across continuous axes using a parametric boundary generator:

| Parameter | Range | Notes |
|---|---|---|
| Boundary family | {linear, additive, bump, enclave, S-shaped, cross, radial} | Plus random combinations |
| # predictive variables | 1–10 | True treatment-effect modifiers |
| # prognostic variables | 0–52 | Affect baseline hazard only |
| Overlap of predictive/prognostic | {none, partial, complete} | |
| Prognostic strength (b₀) | 0–2 | Baseline log-HR |
| Prognostic form | {linear, quadratic, step, saturating} | |
| Censoring rate | 10%–70% | At τ = 30 |
| Covariate correlation | {low, moderate, high} | Block-diagonal or AR(1) |

Target: 200 diverse configurations per n, 5 reps each.

The key insight (from Yue's observation): we don't need many reps per
configuration.  The trend is clear after 3–5 reps.  Diversity of
configurations matters more than precision within a configuration.

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
7. ΔC-index: Cox(Z + A + Z:A) minus Cox(Z)
8. Log-rank test p-value: Cox residuals split by treatment
9. Proportion of β_trt interactions with p < 0.1
10. Interaction F-test p-value from Cox with Z + Z:A

### C. Treatment Effect Profile (bin prognostic score)
11. Variance of RMST difference across 4 bins
12. Slope of RMST difference vs. mean-bin score
13. Deviation from linearity (quadratic contrast) in RMST bins
14. Maximum pairwise RMST diff between bins

### D. Data Quality
15. Censoring rate (overall and per arm)
16. Event count and event rate per arm
17. Mean / median follow-up time
18. Ratio of events to covariates (E/p)

### E. Sample Efficiency Indicators
19. Bootstrap variance of Cox C-index (uncertainty in prognostic signal)
20. Effective sample size (information-based)

### F. Internal Model Behavior (Orig RMSTBoost on training data)
21. Prediction variance (Orig model)
22. StratCF (K=4) prediction variance
23. Pairwise prediction MSE between Orig and StratCF (K=4)
24. Mean subgroup probability (p_i) for each K level
25. Proportion of patients with p_i in ambiguity zone [0.4, 0.6]

### G. Context
26. Sample size n (explicit feature — critical for scaling)
27. Number of covariates p

## Gate Architecture

### Training Labels
For each simulation rep:
1. Fit RMSTBoost for each K ∈ {1, 2, 3, 4, 5} on training data
2. Evaluate AUC on held-out test data (n = 2000)
3. Oracle label = argmax_k AUC_k

### Model
Regularized multinomial logistic regression (LASSO) with softmax output
over K ∈ {1, 2, 3, 4, 5}:

```
P(K = k | x) = exp(β_k · x) / Σ_j exp(β_j · x)
```

where x is the ~40-dimensional feature vector (standardized).

### Optimization
- Loss: categorical cross-entropy
- Regularization: 10-fold CV to select LASSO λ
- The intercept terms encode the *marginal* preference for each K when no
  features are informative (e.g., if the null simulation says K=4 is best
  on average)

### Class Balancing
K=1 and K=5 may be rare in some regimes.  Use balanced accuracy or
weighted loss if needed.

## Validation Strategy

### 1. Within-Family (reproducibility)
80/20 train/test split across all 1,000 runs.  Checks that the gate
learns reproducible patterns.

### 2. Cross-Boundary (generalization)
Train on {linear, additive, bump} families.  Test on {enclave, S-shaped,
cross, radial}.  This is the hardest test: does the gate learn general
data properties, not just boundary-specific lookup tables?

### 3. Null Simulation (no treatment effect)
Gate sees data where A is randomly assigned and has no true effect.
Optimal K should = 1 (no stratification needed when no subgroups exist).
Gate should output K=1 or weight all K equally.

### 4. Comparison Against CV Baseline
For each rep in the test set, compare:
- **Oracle AUC**: best possible (max over K from test data)
- **MoE-K AUC**: gate-chosen K
- **CV AUC**: 5-fold CV chosen K
- **Default AUC**: fixed K=4

The difference (oracle − MoE-K) vs (oracle − CV) is the primary metric.

### 5. Feature Selection Analysis
Which 5–15 features does LASSO retain?
- If ΔC-index is selected → the gate is using interaction strength
- If n is selected → sample size matters for K choice
- If score skewness is selected → enclave detection is happening

This becomes the interpretability result in the paper.

## Data Storage

Each simulation rep is saved as a full RDS file (~0.9–1.2 MB depending
on sample size).  For 1,000 reps, this is approximately 1.2 GB total.

Files stored at `moe/raw/<family>_<seed>.rds` contain:
- `train`: training data (covariates + trt + time + status)
- `test`: test data (covariates + trt + time + status)
- `oracle_label`: TRUE/FALSE for each test patient (treatment benefit)
- `true_te`: true treatment effect on log-hazard scale
- `metadata`: all scenario parameters

This allows recomputing features, fitting new models, or running
additional analyses without re-running the simulation generator.

## Expected Deliverables
1. Parametric scenario sampler (R) ✓
2. Gate feature extraction pipeline (R) — computes ~40 summary stats from training data
3. Full simulation loop (R) — generates data, fits K=1..5, extracts features, records oracle AUC
4. Trained LASSO gate (coefficient vector, ~5–15 non-zero)
5. Validation report:
   - AUC(adaptive K) vs AUC(fixed K=4) vs AUC(CV) vs AUC(oracle)
   - Performance stratified by n
   - Cross-boundary generalization results
   - Null simulation results
   - Selected features and their interpretation

## File Structure
```
moe/
├── README.md                  # This plan
├── scenario_generator.R       # Parametric scenario sampler
├── moe_simulation.R           # Full simulation (data gen → fit → eval → feature extract)
├── gate_features.R            # Feature computation from training data
├── train_gate.R               # LASSO gate training
├── evaluate_gate.R            # Validation on held-out scenarios
├── gate_coefficients.csv      # Trained gate output
└── report.md                  # Validation report
```

## Open Questions
1. How much does adaptive K gain over fixed K=4? (Primary result.)
2. When does MoE-K beat CV? At n=200? n=500? (Secondary result.)
3. Does the gate generalize to boundary families it never saw? (Robustness.)
4. Should K be an integer (discrete choice) or a real-valued K for a
   smoother bias-variance tradeoff? (Design decision.)
