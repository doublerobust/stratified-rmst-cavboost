#!/usr/bin/env Rscript
# gate_importance_viz.R
# Load gate training data, train a ranger RF, extract impurity importance,
# and produce heatmap visualizations showing which feature domains drive gate decisions.
#
# Usage: Rscript moe/gate_importance_viz.R

suppressPackageStartupMessages({
  library(ranger)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(RColorBrewer)
  library(viridisLite)
})

set.seed(20260530)

# ---- Paths ----
DATA_PATH <- "moe/results/gate_training_data.csv"
OUT_DIR   <- "moe/results"

# ---- Feature Domain Groupings (hard-coded) ----
DOMAIN_MAP <- list(
  prognostic_signal       = c("c_index", "score_skew", "score_kurt", "score_var",
                              "score_q90_q10", "orig_pred_var", "orig_pred_mean", "orig_ambiguity"),
  te_heterogeneity        = c("te_bin_var", "te_slope", "te_quadratic_p", "te_max_diff",
                              "te_bin_event_rate_range", "te_int_max_z", "te_int_mean_z", "te_int_prop_sig"),
  interaction_structure   = c("delta_c_index", "prop_interact_sig"),
  data_maturity           = c("censoring_rate", "event_rate", "event_count_trt", "event_count_ctrl",
                              "median_followup", "mean_followup"),
  sample_regime           = c("n", "p", "e_per_p", "cfg_n_predictive", "cfg_n_prognostic",
                              "cfg_te_scale", "cfg_b0"),
  within_arm_cindex       = c("c_index_trt", "c_index_ctrl", "c_index_ratio"),
  correlation             = c("corr_mean", "corr_max", "corr_prop_high"),
  bootstrap_uncertainty   = c("bootstrap_ci_sd", "bootstrap_ci_mean")
)

# Build reverse lookup: feature -> domain
FEATURE_TO_DOMAIN <- list()
for (dm in names(DOMAIN_MAP)) {
  for (ft in DOMAIN_MAP[[dm]]) {
    FEATURE_TO_DOMAIN[[ft]] <- dm
  }
}

# Order of features by domain groups (for consistent column ordering)
DOMAIN_ORDER <- names(DOMAIN_MAP)
ALL_FEATURES <- unlist(DOMAIN_MAP[DOMAIN_ORDER])

# ---- Domain colors ----
DOMAIN_COLORS <- c(
  prognostic_signal     = "#E41A1C",
  te_heterogeneity      = "#377EB8",
  interaction_structure = "#4DAF4A",
  data_maturity         = "#984EA3",
  sample_regime         = "#FF7F00",
  within_arm_cindex     = "#FFFF33",
  correlation           = "#A65628",
  bootstrap_uncertainty = "#F781BF"
)

# ---- Load data ----
cat("Reading", DATA_PATH, "...\n")
df <- read.csv(DATA_PATH, stringsAsFactors = FALSE)
cat(sprintf("Loaded: %d rows x %d cols\n", nrow(df), ncol(df)))

# Columns to exclude from features
exclude_cols <- c("seed", "family", "n_train", "cfg_prognostic_form",
                  "optimal_method", "auc_K1", "auc_K2", "auc_K3", "auc_K4", "auc_K5", "auc_VT")
exclude_cols <- intersect(exclude_cols, names(df))

feature_cols <- setdiff(names(df), exclude_cols)
# Only keep columns that are in our domain map
feature_cols <- intersect(feature_cols, ALL_FEATURES)
cat(sprintf("Using %d feature columns\n", length(feature_cols)))

# Check which features from the domain map are missing from data
missing_feats <- setdiff(ALL_FEATURES, feature_cols)
if (length(missing_feats) > 0) {
  cat("Warning: features not found in data:", paste(missing_feats, collapse = ", "), "\n")
}

# ---- Prepare training data ----
X <- df[, feature_cols, drop = FALSE]

# Impute NAs with column means
for (j in seq_len(ncol(X))) {
  col_mean <- mean(X[[j]], na.rm = TRUE)
  if (is.na(col_mean) || is.nan(col_mean)) col_mean <- 0
  X[[j]][is.na(X[[j]])] <- col_mean
  # Replace non-finite with 0
  X[[j]][!is.finite(X[[j]])] <- 0
}

# Target variable
if ("optimal_method" %in% names(df)) {
  y <- as.factor(df$optimal_method)
  cat("Target: optimal_method (factor levels =", paste(levels(y), collapse = ", "), ")\n")
} else {
  stop("No target column found: expected 'optimal_method'")
}

stopifnot(nrow(X) == length(y))

# ---- Step 1: Train ranger RF ----
n_features <- ncol(X)
mtry_val <- floor(sqrt(n_features))
cat(sprintf("Training ranger RF: %d trees, mtry=%d, %d features, %d rows\n",
            500, mtry_val, n_features, nrow(X)))

rf <- ranger(
  dependent.variable.name = "y",
  data = cbind(X, y = y),
  importance = "impurity",
  num.trees = 500,
  mtry = mtry_val,
  probability = TRUE,
  seed = 20260530,
  verbose = TRUE
)

cat("RF OOB prediction error:", rf$prediction.error, "\n")

# ---- Extract importance ----
imp <- sort(rf$variable.importance, decreasing = TRUE)
imp_df <- data.frame(
  feature    = names(imp),
  importance = as.numeric(imp),
  domain     = sapply(names(imp), function(f) {
    if (!is.null(FEATURE_TO_DOMAIN[[f]])) FEATURE_TO_DOMAIN[[f]] else "other"
  }),
  stringsAsFactors = FALSE
)
imp_df$feature <- factor(imp_df$feature, levels = rev(names(imp)))

mean_imp <- mean(imp_df$importance)
cat(sprintf("Mean importance: %.6f\n", mean_imp))

# ---- Step 2: Global Importance Bar Chart (Figure 1) ----
cat("Creating Figure 1: Global importance bar chart ...\n")

p1 <- ggplot(imp_df, aes(x = importance, y = feature, fill = domain)) +
  geom_col(width = 0.85) +
  geom_vline(xintercept = mean_imp, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  scale_fill_manual(
    values = DOMAIN_COLORS,
    name   = "Feature Domain",
    breaks = names(DOMAIN_COLORS)
  ) +
  labs(
    title = "Gate RF: Global Feature Importance by Domain",
    x = "Impurity Importance",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title          = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position     = "right",
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank(),
    axis.text.y         = element_text(size = 9)
  )

ggsave(file.path(OUT_DIR, "gate_global_importance.pdf"),
       p1, width = 8, height = 10, device = "pdf")
cat("  -> Saved gate_global_importance.pdf\n")

# ---- Prepare for activation heatmaps (Step 3 & 4) ----
# Split data into train/test (80/20 stratified by y)
set.seed(20260530)
n <- nrow(X)
train_idx <- sample(seq_len(n), size = floor(0.8 * n))
test_idx  <- setdiff(seq_len(n), train_idx)

X_train <- X[train_idx, , drop = FALSE]
X_test  <- X[test_idx, , drop = FALSE]
y_test  <- y[test_idx]

# Compute train medians and IQRs for each feature
train_median <- sapply(X_train, median, na.rm = TRUE)
train_iqr    <- sapply(X_train, IQR, na.rm = TRUE)
train_iqr[train_iqr == 0] <- 1  # avoid division by zero

# Global importance sqrt
global_imp_sqrt <- sqrt(imp[feature_cols])
global_imp_sqrt[is.na(global_imp_sqrt) | !is.finite(global_imp_sqrt)] <- 0

# Compute activation for each test row
compute_activation <- function(row_vals, med, iqr, imp_sqrt) {
  z <- (row_vals - med) / iqr
  z <- pmax(pmin(z, 3), -3)
  act <- z * imp_sqrt
  names(act) <- names(row_vals)
  act
}

# Build activation matrix: rows = test datasets, cols = features
act_list <- list()
for (i in seq_len(nrow(X_test))) {
  act_list[[i]] <- compute_activation(
    unlist(X_test[i, feature_cols, drop = TRUE]),
    train_median[feature_cols],
    train_iqr[feature_cols],
    global_imp_sqrt[feature_cols]
  )
}
act_mat <- do.call(rbind, act_list)
rownames(act_mat) <- seq_len(nrow(act_mat))

# ---- Step 3: Multi-Dataset Activation Heatmap (Figure 2) ----
cat("Creating Figure 2: Multi-dataset activation heatmap ...\n")

# Test dataset labels: family_n_train
family_vals   <- df[test_idx, "family", drop = TRUE]
n_train_vals  <- df[test_idx, "n_train", drop = TRUE]
dataset_labels <- paste0(family_vals, "_", n_train_vals)

# Order columns by domain groups
col_order <- intersect(ALL_FEATURES, colnames(act_mat))
act_mat <- act_mat[, col_order, drop = FALSE]

# Build annotation: domain for each column
col_domain <- sapply(col_order, function(f) {
  if (!is.null(FEATURE_TO_DOMAIN[[f]])) FEATURE_TO_DOMAIN[[f]] else "other"
})
names(col_domain) <- col_order

# Build long-format data frame
heatmap_rows <- list()
for (i in seq_len(nrow(act_mat))) {
  for (j in seq_len(ncol(act_mat))) {
    heatmap_rows[[length(heatmap_rows) + 1]] <- data.frame(
      dataset    = dataset_labels[i],
      feature    = col_order[j],
      activation = act_mat[i, j],
      domain     = col_domain[j],
      stringsAsFactors = FALSE
    )
  }
}
heat_df <- do.call(rbind, heatmap_rows)
heat_df$dataset <- factor(heat_df$dataset, levels = rev(unique(dataset_labels)))
heat_df$feature <- factor(heat_df$feature, levels = col_order)

# Domain annotation header
domain_annot <- data.frame(
  feature = col_order,
  domain  = col_domain,
  stringsAsFactors = FALSE
)
domain_annot$feature <- factor(domain_annot$feature, levels = col_order)

# Domain separator positions
domain_tbl <- table(col_domain)[unique(col_domain)]
domain_seps <- data.frame(x = cumsum(domain_tbl) + 0.5)
domain_seps <- domain_seps[-nrow(domain_seps), , drop = FALSE]

# Build the heatmap
p2 <- ggplot(heat_df, aes(x = feature, y = dataset, fill = activation)) +
  geom_tile() +
  scale_fill_gradient2(
    low  = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, name = "Activation"
  ) +
  geom_vline(data = domain_seps, aes(xintercept = x),
             color = "grey30", linewidth = 0.4, linetype = "solid") +
  labs(
    title = "Gate RF: Per-Dataset Activation Map (Test Set)",
    x = NULL, y = "Dataset"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y       = element_text(size = 6),
    panel.grid        = element_blank(),
    legend.position   = "right",
    legend.key.height = unit(1.5, "cm")
  )

# Domain colored top bar
domain_bar <- ggplot(domain_annot, aes(x = feature, y = 1, fill = domain)) +
  geom_tile() +
  scale_fill_manual(
    values = DOMAIN_COLORS,
    name   = "Domain",
    breaks = names(DOMAIN_COLORS)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin     = margin(0, 0, 0, 0),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank()
  ) +
  geom_vline(data = domain_seps, aes(xintercept = x),
             color = "grey30", linewidth = 0.4, linetype = "solid")

# Combine domain bar + heatmap
p2_combined <- grid.arrange(
  domain_bar, p2,
  nrow = 2, ncol = 1,
  heights = c(0.3, 4)
)

ggsave(file.path(OUT_DIR, "gate_activation_heatmap.pdf"),
       p2_combined, width = 12, height = 8, device = "pdf")
cat("  -> Saved gate_activation_heatmap.pdf\n")

# ---- Step 4: Single-Dataset "Brain Slice" (Figure 3) ----
cat("Creating Figure 3: Single-dataset brain slice ...\n")

# Find representative test dataset (closest to median feature vector)
test_med <- apply(act_mat, 2, median, na.rm = TRUE)
dists <- apply(act_mat, 1, function(row) sqrt(mean((row - test_med)^2, na.rm = TRUE)))
rep_idx <- which.min(dists)
rep_label <- dataset_labels[rep_idx]
rep_act <- act_mat[rep_idx, ]

cat(sprintf("  Representative dataset: %s (test row %d)\n", rep_label, rep_idx))

# Domain-level aggregates
domains_in_data <- unique(col_domain)
domain_agg <- data.frame(
  domain = domains_in_data,
  mean_activation = sapply(domains_in_data, function(d) mean(rep_act[col_domain == d])),
  mean_abs_act    = sapply(domains_in_data, function(d) mean(abs(rep_act[col_domain == d]))),
  sign            = sapply(domains_in_data, function(d) {
    val <- mean(rep_act[col_domain == d])
    ifelse(val > 0, "Positive", ifelse(val < 0, "Negative", "Zero"))
  }),
  stringsAsFactors = FALSE
)
domain_agg$domain <- factor(domain_agg$domain, levels = DOMAIN_ORDER)

# Feature-level 1-row strip
strip_df <- data.frame(
  feature    = col_order,
  activation = rep_act,
  domain     = col_domain,
  stringsAsFactors = FALSE
)
strip_df$feature <- factor(strip_df$feature, levels = col_order)

# Domain header positions
domain_tbl_raw <- table(col_domain)[unique(col_domain)]
domain_tbl_features <- as.numeric(domain_tbl_raw)
domain_names         <- names(domain_tbl_raw)
domain_cum <- cumsum(domain_tbl_features)
domain_mid <- c(1, head(domain_cum, -1) + 1) + domain_tbl_features / 2 - 0.5
domain_cum_seps <- head(domain_cum, -1) + 0.5

# Build the brain slice plot
p3_top <- ggplot(strip_df, aes(x = feature, y = 1, fill = activation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, name = "Activation"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = bquote("Gate RF: Activation Map \u2014" ~ .(rep_label))) +
  theme_void() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "none",
    plot.margin     = margin(15, 5, 0, 5),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank()
  ) +
  annotate("text", x = domain_mid, y = 1.4,
           label = gsub("_", " ", domain_names),
           size = 2.8, angle = 0, hjust = 0.5, fontface = "bold") +
  geom_vline(
    data = data.frame(x = domain_cum_seps),
    aes(xintercept = x), color = "grey40", linewidth = 0.5, linetype = "dashed"
  )

# Bottom: domain-level aggregate bars
p3_bottom <- ggplot(domain_agg, aes(x = domain, y = mean_abs_act, fill = sign)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(
    values = c("Positive" = "#B2182B", "Negative" = "#2166AC", "Zero" = "grey80"),
    name   = "Direction"
  ) +
  labs(x = "Domain", y = "Mean |Activation|") +
  theme_minimal(base_size = 9) +
  theme(
    plot.title          = element_blank(),
    axis.text.x         = element_text(angle = 30, hjust = 1, size = 8, face = "bold"),
    axis.text.y         = element_text(size = 7),
    panel.grid.major.x  = element_blank(),
    panel.grid.minor    = element_blank(),
    legend.position     = "right",
    plot.margin         = margin(5, 5, 5, 5)
  )

p3_combined <- grid.arrange(
  p3_top, p3_bottom,
  nrow = 2, ncol = 1,
  heights = c(1, 1.2)
)

ggsave(file.path(OUT_DIR, "gate_brain_slice.pdf"),
       p3_combined, width = 14, height = 4, device = "pdf")
cat("  -> Saved gate_brain_slice.pdf\n")

cat("\nDone! All figures saved to", OUT_DIR, "\n")
