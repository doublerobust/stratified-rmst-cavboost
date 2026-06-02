# Pre-Trained Meta-Learning for Zero-shot Structural Decisions in Statistical Modeling: A Robust Alternative to Cross Validation

## 1. Introduction

### 1.1 The Limitation of Cross-Validation at Small Sample Sizes

Cross-validation is the workhorse of model selection. It estimates out-of-sample performance using within-sample data, and for large samples it works well. But in the small-sample regime that characterizes many clinical trials, genomic studies, and biomedical settings — where n = 200–500 — cross-validation has a fundamental limitation: the standard error of a cross-validated performance estimate is comparable to or larger than the true performance difference between competing methods [Kohavi, 1995; Bengio & Grandvalet, 2004; Varma & Simon, 2006; Hastie et al., 2009, §7.10].

This problem is not restricted to any particular modeling task. It affects hyperparameter tuning (how many trees? which shrinkage rate?), model class selection (boosting vs. random forest vs. neural net?), and structural decisions (should we stratify? what cutpoint?). In all these cases, the practitioner must make a choice based on data that is too sparse to inform it reliably. Cross-validation, applied to a single dataset, cannot see the forest for the trees: it has no access to information from other datasets, other trials, or other studies that faced similar structural decisions.

### 1.2 A Pre-Training Alternative

We investigate a fundamentally different approach. Instead of relying on per-dataset cross-validation, we pre-train a meta-model — a "gate" — on thousands of synthetic datasets that collectively span the range of structural scenarios a practitioner might encounter. This gate learns to map dataset-level summary statistics (the "landmarking" features familiar from the meta-learning literature) to optimal structural decisions. Once pre-trained, the gate applies to a new dataset in a single forward pass: compute a small set of summary statistics (a few seconds of computation), feed them through the frozen gate, and obtain a recommendation.

This is, in spirit, what foundation models do: pre-train on a broad distribution, transfer zero-shot to new instances. The scale in our investigation is smaller — thousands of synthetic datasets rather than billions of tokens — but the principle is the same. The critical difference is that our pre-training is performed on *interpretable meta-features*, not raw covariates. When the gate recommends a particular decision, it does so because the dataset's prognostic strength, interaction signal, and sample size resemble those of other datasets where that decision was optimal.

### 1.3 The Worked Example: Selecting a Subgroup Identification Method

We instantiate this framework with a concrete example: selecting among six competing methods for time-to-event subgroup identification. The candidate methods are:

* **Stratified RMST boosting with K = 1 through 5 strata** (five variants differing in how coarsely or finely the prognostic score is partitioned), and
* **Virtual Twins** (a random-survival-forest-based approach that directly models the treatment-covariate interaction).

Each of the six candidates occupies a different point on the bias-variance spectrum. Low-K stratification (K = 1, i.e., no stratification) is simple and stable but may leave prognostic confounding unaddressed. Higher K values refine the subgroup boundary at the cost of estimating treatment effects within smaller strata. Virtual Twins is the most flexible — modeling interactions nonparametrically — but is also the most data-demanding and can overfit severely at small sample sizes.

A practitioner facing this choice must weigh these tradeoffs using the data at hand. Cross-validation would estimate each candidate's performance on held-out folds, but as argued above, those estimates are too noisy to distinguish among methods when n = 200–500. The gate resolves this by leveraging meta-features — dataset-level summary statistics computed from the training data alone — to predict which candidate will perform best.

Our results show that the pre-trained gate outperforms per-dataset 5-fold cross-validation across all sample sizes. At n = 200, the gate achieves a mean AUC of 0.670 versus 0.655 for CV — a 1.5-point improvement. At n = 300, the gap is 1.6 points (0.698 vs. 0.682). At n = 500 and n = 1000, the gap narrows to 0.7 points (0.773 vs. 0.767 and 0.846 vs. 0.840, respectively). In head-to-head comparisons, the gate beats CV in 28.8% of scenarios versus 19.9% for CV, with 51.4% ties — a 9-percentage-point advantage. The gate uses 26 interpretable dataset-level features, with the prognostic C-index, the bootstrap stability of that estimate, and the quadratic deviation of the treatment-effect profile emerging as the most informative predictors.

### 1.4 The Larger Thesis

The stratified RMST example is a demonstration, not the destination. The meta-learning framework we describe — (a) identify a structural decision that CV handles poorly at small n, (b) generate a broad simulation corpus that covers the decision's input space, (c) pre-train a gate on interpretable dataset-level features, (d) deploy zero-shot — applies across statistical modeling. Any setting where cross-validation for model selection is unreliable at realistic sample sizes is a candidate.

We are not proposing a new method for stratification, nor a new approach to subgroup identification. We are proposing a new paradigm for how structural decisions can be made in data-poor environments: by leveraging information from synthetic data that mimics the breadth of real-world scenarios, summarized into interpretable meta-features that a gate can learn from.

### 1.5 Contributions

1. **A meta-learning framework for structural decisions**: We formalize the problem of replacing per-dataset cross-validation with a pre-trained gate operating on dataset-level meta-features.

2. **A worked example via method selection for subgroup identification**: We demonstrate the framework on the concrete problem of selecting among six candidate methods (K = 1–5 stratified RMST boosting and Virtual Twins), showing that the gate outperforms cross-validation at small n.

3. **A landmarking feature library for survival data**: Twenty-six interpretable features covering prognostic signal, interaction structure, treatment-effect profile, and data quality.

4. **Empirical evidence for the small-n regime**: The gate beats CV in 28.8% of head-to-head comparisons versus 19.9% for CV (51.4% ties), with consistent advantages across n = 200–1000.

5. **Interpretability**: The random forest gate uses impurity importance across all features, providing insight into which data properties drive optimal structural decisions.

### 1.6 Organization

Section 2 reviews related work on cross-validation in small samples, meta-learning for algorithm selection, and the connection to foundation model pre-training. Section 3 describes the meta-learning pipeline: simulation generation, landmarking features, and gate architecture. Section 4 presents results on the stratified RMST example, including the comparison with cross-validation. Section 5 discusses which features the gate learns and the broader applicability of the framework. Section 6 concludes.

## References

1. Kohavi, R. (1995). A study of cross-validation and bootstrap for accuracy estimation and model selection. *Proceedings of the 14th International Joint Conference on Artificial Intelligence* (IJCAI), 2, 1137–1143.

2. Bengio, Y. & Grandvalet, Y. (2004). No unbiased estimator of the variance of K-fold cross-validation. *Journal of Machine Learning Research*, 5, 1089–1105.

3. Varma, S. & Simon, R. (2006). Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*, 7, 91. doi:10.1186/1471-2105-7-91

4. Hastie, T., Tibshirani, R. & Friedman, J. (2009). *The Elements of Statistical Learning* (2nd ed.). Springer. doi:10.1007/978-0-387-84858-7
