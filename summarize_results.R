# Summarize simulation results
# Reads CSV files from results_csv/ and computes mean AUC by scenario

.libPaths(c(file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6"), .libPaths()))

library(dplyr)

csv_dir <- "results_csv"
files <- list.files(csv_dir, pattern = "Sc[0-9]+_Rep[0-9]+\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No CSV result files found in ", csv_dir)
}

cat("Reading", length(files), "result files...\n")

results <- bind_rows(lapply(files, read.csv))

cat(sprintf("Scenarios: %d  (1-%d)\n", length(unique(results$sc)), max(results$sc)))
cat(sprintf("Reps: %d  (range: %d-%d)\n", length(unique(results$rep)), min(results$rep), max(results$rep)))
cat(sprintf("Total rows: %d\n\n", nrow(results)))

# Summary by scenario
summary <- results %>%
  mutate(sc = as.character(sc)) %>%
  group_by(sc) %>%
  summarise(
    n = n(),
    orig_mean = mean(orig_auc, na.rm = TRUE),
    orig_se   = sd(orig_auc, na.rm = TRUE) / sqrt(n()),
    strat_mean = mean(strat_auc, na.rm = TRUE),
    strat_se   = sd(strat_auc, na.rm = TRUE) / sqrt(n()),
    vt_mean   = mean(vt_auc, na.rm = TRUE),
    vt_se     = sd(vt_auc, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Overall average
overall <- results %>%
  summarise(
    n = n(),
    orig_mean = mean(orig_auc, na.rm = TRUE),
    orig_se   = sd(orig_auc, na.rm = TRUE) / sqrt(n()),
    strat_mean = mean(strat_auc, na.rm = TRUE),
    strat_se   = sd(strat_auc, na.rm = TRUE) / sqrt(n()),
    vt_mean   = mean(vt_auc, na.rm = TRUE),
    vt_se     = sd(vt_auc, na.rm = TRUE) / sqrt(n())
  )
overall$sc <- "All"

summary <- bind_rows(summary, overall)

cat("=== Mean AUC by scenario ===\n")
cat(sprintf("%-5s  %-5s  %-20s  %-20s  %-20s\n",
            "Sc", "N", "Orig (SE)", "StratCF (SE)", "VT (SE)"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (i in seq_len(nrow(summary))) {
  r <- summary[i, ]
  cat(sprintf("%-5s  %-5d  %6.4f (%5.4f)  %6.4f (%5.4f)  %6.4f (%5.4f)\n",
              r$sc, r$n,
              r$orig_mean, r$orig_se,
              r$strat_mean, r$strat_se,
              r$vt_mean, r$vt_se))
}

# Pairwise win counts: how often each method has the highest AUC
cat("\n=== Win counts (times each method has highest AUC) ===\n")
wins <- results %>%
  rowwise() %>%
  mutate(
    best = which.max(c(orig_auc, strat_auc, vt_auc)),
    best_method = c("Orig", "StratCF", "VT")[best]
  ) %>%
  ungroup() %>%
  count(best_method) %>%
  rename(Method = best_method, Wins = n)
print(as.data.frame(wins), row.names = FALSE)

cat(sprintf("\nOutput directory: %s/\n", csv_dir))