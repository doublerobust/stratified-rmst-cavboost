#' Gate features: extract dataset-level summary statistics from training data
#' 
#' All features are computed in one pass on the training data — no CV needed.
#' These ~40 features feed into the LASSO gate that predicts optimal K.

library(survival)
library(glmnet)

# =========================================================================
# A. Prognostic Signal (fit Cox on Z only, no treatment)
# =========================================================================

.extract_prognostic_features <- function(data, tau) {
  # Identify covariate columns (exclude trt01p, time, status)
  zcols <- setdiff(names(data), c("trt01p", "time", "status"))
  X <- as.matrix(data[, zcols, drop = FALSE])
  
  # Fit Cox model on covariates only (no treatment)
  cox_null <- tryCatch(
    coxph(Surv(time, status) ~ ., data = data[, c("time", "status", zcols)], 
          iter.max = 25),
    error = function(e) NULL
  )
  
  if (is.null(cox_null)) {
    return(c(
      c_index = NA, score_skew = NA, score_kurt = NA,
      score_var = NA, score_q90_q10 = NA
    ))
  }
  
  # C-index
  c_idx <- summary(cox_null)$concordance[1]
  
  # Prognostic score (linear predictor)
  lp <- predict(cox_null, newdata = data[, zcols], type = "lp")
  lp_centered <- lp - mean(lp)
  
  c(
    c_index = as.numeric(c_idx),
    score_skew = as.numeric(moments::skewness(lp_centered)),
    score_kurt = as.numeric(moments::kurtosis(lp_centered)),
    score_var = as.numeric(var(lp_centered)),
    score_q90_q10 = as.numeric(quantile(lp, 0.9) / quantile(lp, 0.1))
  )
}

# =========================================================================
# B. Interaction Structure (Cox with vs without treatment interactions)
# =========================================================================

.extract_interaction_features <- function(data) {
  zcols <- setdiff(names(data), c("trt01p", "time", "status"))
  
  # Cox: time ~ Z (no treatment) — for baseline C-index
  cox_null <- tryCatch(
    coxph(Surv(time, status) ~ ., data = data[, c("time", "status", zcols)], 
          iter.max = 25),
    error = function(e) NULL
  )
  
  # Main-effects model with treatment (for delta_c_index and trt_main_p)
  main_formula <- as.formula(
    paste("Surv(time, status) ~ trt01p +", paste(zcols, collapse = "+"))
  )
  cox_main <- tryCatch(
    coxph(main_formula, data = data, iter.max = 25),
    error = function(e) NULL
  )
  
  if (is.null(cox_null) || is.null(cox_main)) {
    return(c(delta_c_index = NA, prop_interact_sig = NA, trt_main_p = NA))
  }
  
  # C-index difference (main effects model minus null)
  c_null <- as.numeric(summary(cox_null)$concordance[1])
  c_main <- as.numeric(summary(cox_main)$concordance[1])
  delta_ci <- c_main - c_null
  
  # Treatment main effect p-value
  trt_sum <- coef(summary(cox_main))
  trt_idx <- which(rownames(trt_sum) == "trt01p")
  trt_p <- if (length(trt_idx) == 1) trt_sum[trt_idx, "Pr(>|z|)"] else NA
  
  # Univariate interaction screening: test A:Z_j one at a time
  sig_count <- 0
  total_tested <- 0
  for (zc in zcols) {
    uni_formula <- as.formula(paste("Surv(time, status) ~", zc, "* trt01p"))
    uni_fit <- tryCatch(
      coxph(uni_formula, data = data, iter.max = 25),
      error = function(e) NULL
    )
    if (is.null(uni_fit)) next
    uni_sum <- coef(summary(uni_fit))
    interact_row <- grep(paste0("^", zc, ":trt01p$|^trt01p:", zc, "$"), rownames(uni_sum))
    if (length(interact_row) == 1) {
      total_tested <- total_tested + 1
      if (uni_sum[interact_row, "Pr(>|z|)"] < 0.1) sig_count <- sig_count + 1
    }
  }
  prop_sig <- if (total_tested > 0) sig_count / total_tested else NA
  
  c(delta_c_index = delta_ci,
    prop_interact_sig = prop_sig,
    trt_main_p = trt_p)
}

# =========================================================================
# C. Treatment Effect Profile (binned by prognostic score)
# =========================================================================

.extract_te_profile_features <- function(data, tau, n_bins = 4) {
  zcols <- setdiff(names(data), c("trt01p", "time", "status"))
  A <- data$trt01p
  
  # Get prognostic score from null Cox
  cox_null <- tryCatch(
    coxph(Surv(time, status) ~ ., data = data[, c("time", "status", zcols)], 
          iter.max = 25),
    error = function(e) NULL
  )
  
  if (is.null(cox_null)) {
    return(c(te_bin_var = NA, te_slope = NA, te_quadratic = NA,
             te_max_diff = NA, te_bin_event_rate_range = NA))
  }
  
  lp <- predict(cox_null, newdata = data[, zcols], type = "lp")
  bins <- cut(lp, breaks = unique(quantile(lp, probs = seq(0, 1, 1/n_bins), na.rm = TRUE)),
              include.lowest = TRUE)
  
  # Compute RMST difference per bin
  rmst_diff <- numeric(n_bins)
  event_rate <- numeric(n_bins)
  bin_center <- numeric(n_bins)
  
  for (k in seq_len(n_bins)) {
    idx <- which(as.integer(bins) == k)
    if (length(idx) < 10) next
    
    bin_center[k] <- mean(lp[idx], na.rm = TRUE)
    event_rate[k] <- mean(data$status[idx])
    
    # RMST in treatment vs control within this bin
    for (arm in 0:1) {
      arm_idx <- idx[which(A[idx] == arm)]
      if (length(arm_idx) < 5) next
      km <- survfit(Surv(time, status) ~ 1, data = data[arm_idx, ])
      # RMST at tau
      sfit <- summary(km, times = tau, extend = TRUE)
      if (arm == 0) rmst0 <- sfit$table["median"] else rmst1 <- sfit$table["median"]
      # Actually compute properly
    }
    # Simple Kaplan-Meier RMST approximation
    km0 <- survfit(Surv(time, status) ~ 1, data = data[idx[A[idx] == 0], , drop = FALSE])
    km1 <- survfit(Surv(time, status) ~ 1, data = data[idx[A[idx] == 1], , drop = FALSE])
    
    # Area under KM curve up to tau
    sfit0 <- summary(km0, times = seq(0, tau, length.out = 100))
    sfit1 <- summary(km1, times = seq(0, tau, length.out = 100))
    
    # Integrate using trapezoidal rule
    rmst0 <- sum(diff(sfit0$time) * (head(sfit0$surv, -1) + tail(sfit0$surv, -1)) / 2)
    rmst1 <- sum(diff(sfit1$time) * (head(sfit1$surv, -1) + tail(sfit1$surv, -1)) / 2)
    
    rmst_diff[k] <- rmst1 - rmst0
  }
  
  valid <- is.finite(rmst_diff) & is.finite(bin_center) & bin_center != 0
  
  if (sum(valid) < 3) {
    return(c(te_bin_var = NA, te_slope = NA, te_quadratic = NA,
             te_max_diff = NA, te_bin_event_rate_range = NA))
  }
  
  # Variance of RMST difference across bins
  te_var <- var(rmst_diff[valid], na.rm = TRUE)
  
  # Linear slope
  te_lm <- lm(rmst_diff[valid] ~ bin_center[valid])
  te_slope_coef <- coef(te_lm)[2]
  
  # Quadratic deviation
  te_lm_quad <- lm(rmst_diff[valid] ~ poly(bin_center[valid], 2))
  te_anova <- anova(te_lm, te_lm_quad)
  te_quad_p <- te_anova$"Pr(>F)"[2]
  
  # Max pairwise difference
  te_max_abs_diff <- if (sum(valid) >= 2) max(dist(rmst_diff[valid])) else NA
  
  # Event rate range across bins
  event_range <- diff(range(event_rate[event_rate > 0], na.rm = TRUE))
  
  c(te_bin_var = as.numeric(te_var),
    te_slope = as.numeric(te_slope_coef),
    te_quadratic_p = as.numeric(te_quad_p),
    te_max_diff = as.numeric(te_max_abs_diff),
    te_bin_event_rate_range = as.numeric(event_range))
}

# =========================================================================
# D. Data Quality
# =========================================================================

.extract_quality_features <- function(data, tau) {
  A <- data$trt01p
  zcols <- setdiff(names(data), c("trt01p", "time", "status"))
  p <- length(zcols)
  
  c(
    censoring_rate = 1 - mean(data$status),
    event_rate = mean(data$status),
    event_count_trt = sum(data$status[A == 1]),
    event_count_ctrl = sum(data$status[A == 0]),
    median_followup = median(data$time, na.rm = TRUE),
    mean_followup = mean(data$time, na.rm = TRUE),
    n = nrow(data),
    p = p,
    e_per_p = sum(data$status) / p
  )
}

# =========================================================================
# E. Sample Efficiency Indicators
# =========================================================================

.extract_efficiency_features <- function(data) {
  zcols <- setdiff(names(data), c("trt01p", "time", "status"))
  
  # Bootstrap variance of C-index (n_boot = 50 for speed)
  n_boot <- 50
  ci_boot <- numeric(n_boot)
  
  for (b in seq_len(n_boot)) {
    idx <- sample(nrow(data), replace = TRUE)
    boot_data <- data[idx, ]
    
    cox_boot <- tryCatch(
      coxph(Surv(time, status) ~ ., data = boot_data[, c("time", "status", zcols)],
            iter.max = 25),
      error = function(e) NULL
    )
    ci_boot[b] <- if (!is.null(cox_boot)) summary(cox_boot)$concordance[1] else NA
  }
  
  c(bootstrap_ci_sd = sd(ci_boot, na.rm = TRUE),
    bootstrap_ci_mean = mean(ci_boot, na.rm = TRUE))
}

# =========================================================================
# F. Internal Model Behavior (from Orig RMSTBoost on training data)
# =========================================================================

# These require the model to be fit first — called from the simulation loop
# with the model object and training predictions

# =========================================================================
# G. Context
# =========================================================================



# =========================================================================
# H. Within-Arm Features (separate treatment and control)
# =========================================================================

.extract_within_arm_features <- function(data) {
  zcols <- setdiff(names(data), c("trt01p", "time", "status"))
  A <- data$trt01p
  
  # C-index within each arm
  c_trt <- NA; c_ctrl <- NA
  for (arm in 1:0) {
    sub <- data[A == arm, c("time", "status", zcols)]
    if (nrow(sub) < 20) next
    cox <- tryCatch(
      coxph(Surv(time, status) ~ ., data = sub, iter.max = 25),
      error = function(e) NULL
    )
    if (!is.null(cox)) {
      if (arm == 1) c_trt <- as.numeric(summary(cox)$concordance[1])
      else c_ctrl <- as.numeric(summary(cox)$concordance[1])
    }
  }
  
  c(c_index_trt = c_trt,
    c_index_ctrl = c_ctrl,
    c_index_ratio = c_trt / c_ctrl)
}

# =========================================================================
# I. Univariate Treatment Interaction Signals
# =========================================================================

.extract_te_interaction_features <- function(data) {
  zcols <- setdiff(names(data), c("trt01p", "time", "status"))
  A <- data$trt01p
  
  z_scores <- numeric(length(zcols))
  for (j in seq_along(zcols)) {
    # Fit Cox with treatment x covariate interaction
    form <- as.formula(paste0("Surv(time, status) ~ trt01p + ", zcols[j], " + trt01p:", zcols[j]))
    cox_int <- tryCatch(
      coxph(form, data = data, iter.max = 25),
      error = function(e) NULL
    )
    if (!is.null(cox_int)) {
      # Extract interaction coefficient z-score
      int_coef <- coef(cox_int)[paste0("trt01p:", zcols[j])]
      if (!is.na(int_coef)) {
        se <- sqrt(diag(vcov(cox_int)))[paste0("trt01p:", zcols[j])]
        z_scores[j] <- abs(int_coef / se)
      } else {
        z_scores[j] <- 0
      }
    } else {
      z_scores[j] <- 0
    }
  }
  
  c(
    te_int_max_z = max(z_scores, na.rm = TRUE),
    te_int_mean_z = mean(z_scores, na.rm = TRUE),
    te_int_prop_sig = mean(z_scores > 1.96, na.rm = TRUE)
  )
}

# =========================================================================
# J. Feature Correlation Structure
# =========================================================================

.extract_correlation_features <- function(data) {
  zcols <- setdiff(names(data), c("trt01p", "time", "status"))
  X <- as.matrix(data[, zcols, drop = FALSE])
  
  cormat <- tryCatch(cor(X, use = "pairwise.complete.obs"), error = function(e) NULL)
  
  if (is.null(cormat) || ncol(cormat) < 2) {
    return(c(corr_mean = NA, corr_max = NA, corr_prop_high = NA))
  }
  
  tril_vals <- cormat[lower.tri(cormat)]
  
  c(
    corr_mean = mean(abs(tril_vals)),
    corr_max = max(abs(tril_vals)),
    corr_prop_high = mean(abs(tril_vals) > 0.5)
  )
}

# =========================================================================
# Master feature extraction function
# =========================================================================

compute_gate_features <- function(data, tau = 30, fit_orig = NULL, model_preds = NULL) {
  features <- c(
    # A. Prognostic
    .extract_prognostic_features(data, tau),
    
    # B. Interaction structure
    .extract_interaction_features(data),
    
    # C. TE profile
    .extract_te_profile_features(data, tau),
    
    # D. Data quality
    .extract_quality_features(data, tau),
    
    # E. Efficiency
    .extract_efficiency_features(data),
    
    # H. Within-arm
    .extract_within_arm_features(data),
    
    # I. TE interaction
    .extract_te_interaction_features(data),
    
    # J. Correlation
    .extract_correlation_features(data),
    
    # n and p explicitly
    n = nrow(data),
    p = length(setdiff(names(data), c("trt01p", "time", "status")))
  )
  
  # F. Model behavior (if provided)
  if (!is.null(fit_orig) && !is.null(model_preds)) {
    pred_var <- var(model_preds, na.rm = TRUE)
    pred_mean <- mean(model_preds, na.rm = TRUE)
    ambiguity <- mean(model_preds > 0.4 & model_preds < 0.6, na.rm = TRUE)
    
    features <- c(features,
      orig_pred_var = pred_var,
      orig_pred_mean = pred_mean,
      orig_ambiguity = ambiguity
    )
  }
  
  features
}
