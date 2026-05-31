#!/usr/bin/env Rscript
#
# Diagnose AUC=1 cases: find scenario characteristics & check for data leakage
#

RESULTS_DIR <- "moe/results"

result_files <- list.files(RESULTS_DIR, "rep_.*\\.rds$", full.names = TRUE)
cat(sprintf("Found %d result files\n\n", length(result_files)))

# ---- 1. Identify near-perfect AUC cases ----
perfect_list <- list()
total <- 0

for (f in result_files) {
  r <- readRDS(f)
  aucs <- r$aucs
  max_auc <- max(aucs, na.rm = TRUE)
  
  if (is.finite(max_auc) && max_auc > 0.999) {
    perfect_list[[length(perfect_list) + 1]] <- data.frame(
      file = basename(f),
      family = r$config$family,
      n_train = r$config$n_train,
      b0 = r$config$b0,
      n_predictive = r$config$n_predictive,
      n_prognostic = r$config$n_prognostic,
      overlap = r$config$overlap,
      prognostic_form = r$config$prognostic_form,
      censoring_rate = r$config$censoring_rate,
      corr = r$config$corr,
      auc_K1 = aucs[1], auc_K2 = aucs[2], auc_K3 = aucs[3],
      auc_K4 = aucs[4], auc_K5 = aucs[5],
      max_auc = max_auc,
      oracle_rate = r$oracle_rate,
      stringsAsFactors = FALSE
    )
  }
  total <- total + 1
}

if (length(perfect_list) == 0) {
  cat("No near-perfect AUC (> 0.999) cases found in", total, "results.\n")
  quit(save = "no", status = 0)
}

perfect_df <- do.call(rbind, perfect_list)
cat(sprintf("Found %d / %d results with AUC > 0.999 (%.1f%%)\n\n",
            nrow(perfect_df), total, 100 * nrow(perfect_df) / total))

# ---- 2. Distribution by scenario characteristics ----
cat("=== By family ===\n")
print(table(perfect_df$family))
cat("\n=== By n_train ===\n")
print(table(perfect_df$n_train))
cat("\n=== By overlap ===\n")
print(table(perfect_df$overlap))
cat("\n=== By prognostic_form ===\n")
print(table(perfect_df$prognostic_form))
cat("\n=== By n_predictive ===\n")
print(table(perfect_df$n_predictive))
cat("\n=== By b0 range ===\n")
print(table(cut(perfect_df$b0, c(-Inf, 0.5, 1.0, 1.5, 2.0, Inf))))
cat("\n=== By censoring_rate range ===\n")
print(table(cut(perfect_df$censoring_rate, c(-Inf, 0.1, 0.3, 0.5, 0.7, Inf))))
cat("\n=== By corr ===\n")
print(table(perfect_df$corr))
cat("\n=== By n_prognostic ===\n")
print(table(perfect_df$n_prognostic))

# ---- 3. Which K gives AUC=1 most often? ----
cat("\n=== Which K values have AUC=1 ===\n")
for (k in 1:5) {
  col <- paste0("auc_K", k)
  n <- sum(perfect_df[[col]] > 0.999, na.rm = TRUE)
  cat(sprintf("  K=%d: %d / %d near-perfect (%.1f%% of all perfect cases)\n",
              k, n, nrow(perfect_df), 100 * n / nrow(perfect_df)))
}

# ---- 4. Detailed check: is this legit? ----
cat("\n=== Detailed check: first 5 near-perfect cases ===\n\n")

n_check <- min(5, nrow(perfect_df))
for (i in seq_len(n_check)) {
  f <- perfect_df$file[i]
  r <- readRDS(file.path(RESULTS_DIR, f))
  
  cat(sprintf("Case %d: %s\n", i, f))
  cat(sprintf("  Config: family=%s, n_train=%d, b0=%.2f, censoring=%.2f\n",
              r$config$family, r$config$n_train, r$config$b0, 
              r$config$censoring_rate))
  cat(sprintf("  AUCs: K1=%.4f K2=%.4f K3=%.4f K4=%.4f K5=%.4f\n",
              r$aucs[1], r$aucs[2], r$aucs[3], r$aucs[4], r$aucs[5]))
  cat(sprintf("  Oracle rate: %.3f\n", r$oracle_rate))
  cat(sprintf("  Oracle optimal K: %d\n", r$oracle_optimal_K))
  cat("\n")
}

# ---- 5. Check for suspicious patterns ----
cat("=== Checking for data leakage patterns ===\n\n")

# Check if AUC=1 is concentrated in specific families
cat("AUC > 0.999 rate by family:\n")
all_summary <- data.frame()
for (f in result_files) {
  r <- readRDS(f)
  all_summary <- rbind(all_summary, data.frame(
    family = r$config$family,
    max_auc = max(r$aucs, na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
}
tbl <- table(all_summary$family, all_summary$max_auc > 0.999)
tbl_prop <- prop.table(tbl, 1)
colnames(tbl_prop) <- c("AUC<1", "AUC>0.999")
print(round(tbl_prop * 100, 1))
cat("\n")

cat("AUC > 0.999 rate by n_train:\n")
tbl2 <- table(perfect_df$n_train)
tbl2_all <- table(all_summary$n_train)
for (n in names(tbl2_all)) {
  n_pt <- if (n %in% names(tbl2)) tbl2[n] else 0
  n_all <- tbl2_all[n]
  cat(sprintf("  n=%s: %d / %d perfect (%.1f%%)\n", n, n_pt, n_all, 100 * n_pt / n_all))
}

# ---- 6. Most common combination leading to AUC=1 ----
cat("\n=== Top 10 most common (family, b0, n_train, censoring) combos for AUC=1 ===\n")
if (nrow(perfect_df) > 0) {
  combo <- paste(perfect_df$family, round(perfect_df$b0, 1), 
                 perfect_df$n_train, perfect_df$censoring_rate)
  combo_tbl <- sort(table(combo), decreasing = TRUE)
  print(head(combo_tbl, 10))
}

cat("\nDone.\n")