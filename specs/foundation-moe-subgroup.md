# Foundation-MoE for Subgroup Identification (Spec v0.1)

## Motivation

Current stratified CAVBoost uses fixed Cox PH prognostic strata (K quantiles).
On disjoint/non-monotone boundaries (S4 Enclave), the fixed strata misdirect
the gradient because complex treatment-effect pockets span multiple strata.

Foundation models (TabICL, TabPFN, TabR) produce rich feature representations
that could replace the Cox PH score as input to a **learned gating network**,
making strata adaptive to the treatment-effect boundary.

## Architecture

```
Z_i (covariates)
    │
    ▼
Foundation model (TabICL encoder / TabPFN)
    │  produces representation ϕ(Z_i)
    │
    ├─► Prognostic Head:   S_i = MLP(ϕ(Z_i))   [continuous score → K strata]
    ├─► Subgroup Head:     η_i  = TreeBoost(ϕ(Z_i))  [subgroup prediction]
    │
    │  (Alternative: single shared trunk → two heads)
    │
    ▼
V_strat(p) = Σ_k W_k·d_k^(1) − Σ_k (n_k−W_k)·d_k^(2)
    where strata are defined by quantiles of S_i from the Prognostic Head
```

**Key difference from current approach:** The prognostic score S_i is learned
end-to-end via backprop through V_strat, rather than pre-computed via Cox PH.
This lets the gradient signal which baseline features are predictive and
adapts the strata boundaries to the treatment-effect structure.

## Design Decisions

| Decision | Options | Recommendation |
|----------|---------|---------------|
| Foundation model | TabICL, TabR, TabPFN | TabICL simplest (kNN-based, no training) |
| Gate architecture | MLP (1-2 layers) | Small — learnable params ≪ n |
| Strata count K | 2, 4 | 2 for stability, 4 for granularity |
| Training | Alternating (fix gate, update trees) | Avoids end-to-end gradient issues |
| Software | Python for gate, R for gradient | Or single Python implementation |

## Why Not End-to-End PyTorch

The weighted KM gradient is sequential (cumulative over sorted event times).
GPU acceleration does not help. Keep R/Python hybrid to preserve the verified
gradient implementation.

## Implementation Plan

Phase 1 — Gate as separate model:
1. Train foundation model representation for all patients
2. Fit a small MLP on top to predict S_i (prognostic score)
3. Use S_i to define strata (same quantile cut)
4. Run current stratified CAVBoost with these strata
5. Compare to Cox PH strata

Phase 2 — Alternating optimization:
1. Fix MLP gate → train XGBoost subgroup head
2. Fix XGBoost → update MLP gate to maximize V_strat on holdout
3. Iterate 5-10 alternating rounds

## Related Work

- Xu et al. (2023): Deep learning subgroup identification with survival outcomes
- TabICL (Hollmann et al. 2023): In-context learning for tabular data
- Balazadeh et al. (2025): CausalPFN — prior-based treatment effect estimation

## Open Questions

1. Does TabICL's representation improve over Cox PH at n=500?
   (Need: breast cancer / PBC validation dataset with ~500 obs)
2. Can a 1-layer MLP learn better strata than quantile cuts of a Cox LP?
3. Does adaptive stratification fix S4 Enclave (AUC gap of ~1.5 pts)?
4. What's the Type I error cost of optimizing strata boundaries?

---
Status: Spec v0.1 — parked for future Tabular Foundation work
