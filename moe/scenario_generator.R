#' MoE-K scenario generator
#' 
#' Parametric sampler for diverse treatment-effect boundary geometries.
#' Generalizes the 6 fixed Zhang scenarios into parameterized families.
#' Each call generates one dataset with known oracle labels.

library(mvtnorm)
library(survival)

# =========================================================================
# 1. Boundary family implementations
# =========================================================================

# All families return a numeric vector te (log-hazard scale, per patient).
# Positive te = treatment benefit, negative = harm.
# X is the covariate matrix with columns: Z1..Z4, z5..z50, S1, S2

.te_linear <- function(X, n_predictive, ...) {
  if (n_predictive >= 1L) {
    te <- X[, "S1"]
    if (n_predictive >= 2L) {
      gamma <- runif(1, 0.5, 1.5)
      te <- te + gamma * X[, "S2"]
    }
  } else {
    te <- rep(0, nrow(X))
  }
  te
}

.te_additive <- function(X, n_predictive, ...) {
  if (n_predictive >= 2L) {
    X[, "S1"] - X[, "S2"]
  } else {
    rep(0, nrow(X))
  }
}

.te_bump <- function(X, ...) {
  s <- X[, "S1"]
  2 * ifelse(abs(s) < 0.67, exp(-s^2) - 0.4, exp(-s^2) - 0.8)
}

.te_enclave <- function(X, ...) {
  s1 <- X[, "S1"]; s2 <- X[, "S2"]
  2 * ((-1.07 <= s1 & s1 < 1.07) & (-1.07 <= s2 & s2 < 1.07)) - 1
}

.te_s_shaped <- function(X, ...) {
  s <- X[, "S1"]
  2 * ifelse(s >= 0.67 | (-0.67 <= s & s < 0), 1, 0) - 1
}

.te_cross <- function(X, ...) {
  s1 <- X[, "S1"]; s2 <- X[, "S2"]
  2 * ifelse((s1 >= 0 & s2 >= -0.67) | (s1 < 0 & s2 < -0.67), 1, 0) - 1
}

.te_radial <- function(X, ..., r0 = NULL) {
  s1 <- X[, "S1"]; s2 <- X[, "S2"]
  r <- sqrt(s1^2 + s2^2)
  if (is.null(r0)) r0 <- runif(1, 0.5, 2.0)
  2 * (r < r0) - 1
}

.te_random <- function(X, ...) {
  # Voronoi-based random boundary
  n_regions <- sample(3:8, 1)
  centers <- matrix(runif(n_regions * 2, -3, 3), ncol = 2)
  signs <- sample(c(-1, 1), n_regions, replace = TRUE)
  s1 <- X[, "S1"]; s2 <- X[, "S2"]
  
  te <- numeric(nrow(X))
  for (i in seq_len(n_regions)) {
    dist_i <- sqrt((s1 - centers[i, 1])^2 + (s2 - centers[i, 2])^2)
    if (i == 1L) {
      nearest <- rep(1L, nrow(X))
      min_dist <- dist_i
    } else {
      closer <- dist_i < min_dist
      nearest[closer] <- i
      min_dist[closer] <- dist_i[closer]
    }
  }
  signs[nearest] * runif(1, 0.5, 1.5)
}

BOUNDARY_FUNCTIONS <- list(
  linear    = .te_linear,
  additive  = .te_additive,
  bump      = .te_bump,
  enclave   = .te_enclave,
  s_shaped  = .te_s_shaped,
  cross     = .te_cross,
  radial    = .te_radial,
  random    = .te_random
)

# =========================================================================
# 2. Main generator
# =========================================================================

#' Generate one simulation scenario
#'
#' @param n_train Training set size
#' @param n_test  Test set size
#' @param family Boundary family name
#' @param n_predictive Number of true predictive variables (1 or 2 for families
#'   that use S1,S2; 1 for S1-only families)
#' @param n_prognostic Number of prognostic variables (0-52)
#' @param overlap Overlap type: "none", "partial", "complete"
#' @param b0 Baseline log-hazard (prognostic strength)
#' @param prognostic_form "linear", "quadratic", or "mixed"
#' @param censoring_rate Target censoring rate at tau (0-1)
#' @param tau Restricted mean survival horizon
#' @param corr Covariate correlation: "low" (0.1), "moderate" (1/3), "high" (0.7)
#' @param n_vars Total number of covariates (default 52)
#' @param seed Random seed
#' @param save_dir If non-NULL, save the full output as RDS at save_dir/family_seed.rds
#' @param ... Additional arguments passed to boundary functions
#'
#' @return List with train, test, oracle_label, true_te, metadata, or NULL if invalid

# Minimum n_predictive required per family
.NEEDS_2_PREDICTIVE <- c("additive", "enclave", "cross", "radial", "random")
.NEEDS_1_PREDICTIVE <- c("linear", "bump", "s_shaped")

generate_scenario <- function(
    n_train = 500,
    n_test = 2000,
    family = "linear",
    n_predictive = if (family %in% .NEEDS_2_PREDICTIVE) 2L else 1L,
    n_prognostic = 4L,
    overlap = "none",
    b0 = sqrt(6),
    prognostic_form = "linear",
    censoring_rate = 0.10,
    tau = 30,
    corr = "moderate",
    n_vars = 52L,
    te_scale = 1.0,
    seed = NULL,
    save_dir = NULL,
    ...
) {
  if (!is.null(seed)) set.seed(seed)
  
  # ---- Validate ----
  stopifnot(
    family %in% names(BOUNDARY_FUNCTIONS),
    n_predictive >= 1L, n_predictive <= 10L,
    n_prognostic >= 0L, n_prognostic <= n_vars - 2L,
    n_train >= 50L, n_test >= 100L
  )
  
  # Validate n_predictive is sufficient for this family
  if (family %in% .NEEDS_2_PREDICTIVE && n_predictive < 2L) {
    stop(sprintf("Family '%s' requires n_predictive >= 2 (got %d)", family, n_predictive))
  }
  
  n_total <- n_train + n_test
  
  # ---- Covariance matrix ----
  rho <- switch(corr,
    low = 0.1,
    moderate = 1/3,
    high = 0.7,
    1/3
  )
  Smat <- matrix(rho, n_vars, n_vars)
  diag(Smat) <- 1
  
  # ---- Generate covariates ----
  X <- rmvnorm(n_total, sigma = Smat)
  colnames(X) <- c(paste0("z", 1:(n_vars - 2L)), "S1", "S2")
  colnames(X)[1:4] <- c("Z1", "Z2", "Z3", "Z4")
  
  # ---- Treatment assignment ----
  A <- rbinom(n_total, 1, 0.5)
  
  # ---- Treatment effect boundary ----
  te <- BOUNDARY_FUNCTIONS[[family]](X, n_predictive = n_predictive, ...)
  
  # Scale treatment effect to realistic level
  te <- te * te_scale
  
  # ---- Oracle label ----
  oracle_label <- as.logical(te > 0)
  if (length(unique(oracle_label)) < 2L) {
    return(NULL)  # invalid boundary for this seed
  }
  
  # ---- Prognostic effect ----
  # Prognostic variables: when overlap = "none", use Z1..Zn_prognostic
  # When overlap = "complete", use S1,S2 + Z1..Zn_prognostic-2
  # When overlap = "partial", mix
  
  if (overlap == "none") {
    prog_vars <- grep("^Z[0-9]", colnames(X))  # Z1..Z4
    if (n_prognostic < 4) prog_vars <- prog_vars[seq_len(n_prognostic)]
  } else if (overlap == "complete") {
    # Predictive + prognostic use same vars: S1,S2 + some Zs
    prog_vars <- c(grep("^S[0-9]", colnames(X)), grep("^Z[0-9]", colnames(X)))
    n_prog_avail <- length(prog_vars)
    prog_vars <- prog_vars[seq_len(min(n_prognostic, n_prog_avail))]
  } else {
    # "partial": S1,S2 are predictive, Z1..Z3 are additionally prognostic
    prog_vars <- c(grep("^S[0-9]", colnames(X)), grep("^Z[0-9]", colnames(X))[1:3])
    prog_vars <- prog_vars[seq_len(min(n_prognostic + n_predictive, length(prog_vars)))]
  }
  
  # Prognostic coefficients
  beta_raw <- runif(length(prog_vars), 0.2, 0.6)
  beta_signs <- sample(c(-1, 1), length(prog_vars), replace = TRUE)
  beta_prog <- beta_raw * beta_signs
  
  Zb <- as.numeric(X[, prog_vars, drop = FALSE] %*% beta_prog)
  
  if (prognostic_form == "quadratic" || 
      (prognostic_form == "mixed" && runif(1) < 0.5)) {
    # Saturating prognostic effect (as in Zhang S4-S6)
    Zb <- -Zb^2
  }
  
  # ---- Generate survival times ----
  s0 <- 0.4
  T <- exp(b0 + A * te + Zb + s0 * rnorm(n_total))
  
  # Censoring: tune rate to hit target
  if (censoring_rate > 0) {
    # Solve for censoring rate parameter approximately
    C <- pmin(tau, rexp(n_total, rate = -log(1 - censoring_rate) / (tau / 2)))
  } else {
    C <- rep(tau, n_total)
  }
  
  U <- pmin(T, C, tau)
  status <- as.numeric(T <= C & U < tau)
  
  # ---- Split into train/test ----
  train_idx <- seq_len(n_train)
  test_idx <- seq(n_train + 1, n_total)
  
  train <- as.data.frame(X[train_idx, , drop = FALSE])
  train$trt01p <- A[train_idx]
  train$time <- U[train_idx]
  train$status <- status[train_idx]
  
  test <- as.data.frame(X[test_idx, , drop = FALSE])
  test$trt01p <- A[test_idx]
  test$time <- U[test_idx]
  test$status <- status[test_idx]
  
  result <- list(
    train = train,
    test = test,
    oracle_label = oracle_label[test_idx],
    train_oracle_label = oracle_label[train_idx],
    true_te = te[test_idx],
    metadata = list(
      family = family,
      n_predictive = n_predictive,
      n_prognostic = n_prognostic,
      te_scale = te_scale,
      overlap = overlap,
      b0 = b0,
      prognostic_form = prognostic_form,
      censoring_rate = censoring_rate,
      corr = corr,
      n_train = n_train,
      n_test = n_test,
      seed = seed,
      prog_vars = colnames(X)[prog_vars],
      beta_prog = beta_prog
    )
  )
  
  # Save raw data if requested
  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    save_path <- file.path(save_dir, sprintf("%s_%d.rds", family, seed))
    saveRDS(result, save_path)
  }
  
  result
}

# =========================================================================
# 3. Scenario configuration sampler
# =========================================================================

#' Sample random scenario configurations for batch simulation
#'
#' @param n_configs Number of configuration samples
#' @return data.frame with one row per configuration
sample_configurations <- function(n_configs = 200) {
  
  families <- c("linear", "additive", "bump", "enclave", 
                "s_shaped", "cross", "radial", "random")
  
  # Which families need n_predictive = 2?
  needs_2 <- c("additive", "enclave", "cross", "radial", "random")
  needs_1 <- c("linear", "bump", "s_shaped")
  
  sampled_families <- sample(families, n_configs, replace = TRUE)
  
  # Assign n_predictive based on family
  n_predictive <- integer(n_configs)
  n_predictive[sampled_families %in% needs_2] <- 2L
  n_predictive[sampled_families %in% needs_1] <- 1L
  n_predictive[!sampled_families %in% c(needs_2, needs_1)] <- sample(1:2, sum(!sampled_families %in% c(needs_2, needs_1)), replace = TRUE)
  
  data.frame(
    family = sampled_families,
    n_predictive = n_predictive,
    n_prognostic = sample(c(0, 1, 2, 4, 8, 16), n_configs, replace = TRUE,
                          prob = c(0.1, 0.1, 0.2, 0.3, 0.2, 0.1)),
    overlap = sample(c("none", "partial", "complete"), n_configs, replace = TRUE,
                     prob = c(0.2, 0.5, 0.3)),
    te_scale = round(runif(n_configs, 0.25, 0.50), 2),
    b0 = round(runif(n_configs, 0, 2), 2),
    prognostic_form = sample(c("linear", "quadratic", "mixed"), n_configs,
                             replace = TRUE, prob = c(0.4, 0.3, 0.3)),
    censoring_rate = round(runif(n_configs, 0.05, 0.70), 2),
    corr = sample(c("low", "moderate", "high"), n_configs, replace = TRUE,
                  prob = c(0.3, 0.4, 0.3)),
    n_train = sample(c(200, 300, 500, 1000), n_configs, replace = TRUE,
                     prob = c(0.2, 0.2, 0.4, 0.2)),
    seed = sample(1:1e6, n_configs, replace = FALSE),
    stringsAsFactors = FALSE
  )
}
