#!/usr/bin/env Rscript
# Compute classification metrics from raw prediction files.
# Reads results_raw/ (agent) and results_raw_omen/ (Omen).
# Outputs summary table.

library(parallel)
options(mc.cores = 6)

raw_dirs <- c("results_raw", "results_raw_omen")
sc_names <- c("S1 Linear", "S2 Diff", "S3 U-shaped",
              "S4 Enclave", "S5 S-shaped", "S6 Cross")

classif <- function(pred, label, cutoff = 0.5) {
  pred_class <- pred > cutoff
  tp <- sum(pred_class & label); fp <- sum(pred_class & !label)
  tn <- sum(!pred_class & !label); fn <- sum(!pred_class & label)
  c(sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    accuracy = (tp + tn) / (tp + fp + tn + fn))
}

process_file <- function(f) {
  r <- readRDS(f)
  label <- as.logical(r$label)
  n <- sum(label)
  orig <- classif(r$orig_pred, label)
  strat <- classif(r$strat_pred, label)
  vt <- classif(r$vt_pred, label)
  data.frame(sc = r$scenario, rep = r$rep,
             orig_sens = orig[1], orig_spec = orig[2], orig_acc = orig[3],
             strat_sens = strat[1], strat_spec = strat[2], strat_acc = strat[3],
             vt_sens = vt[1], vt_spec = vt[2], vt_acc = vt[3])
}

cat("Reading raw files...\n")
all_files <- unlist(lapply(raw_dirs, function(d) {
  fs <- list.files(d, "Sc.*_raw\\.rds$", full.names = TRUE)
  cat(sprintf("  %s: %d files\n", d, length(fs)))
  fs
}))

res <- do.call(rbind, lapply(all_files, process_file))
cat(sprintf("Total: %d reps\n", nrow(res)))

cat("\n=== Classification Metrics ===\n\n")
cat(sprintf("%-14s %-5s %-10s %-10s %-10s %-10s %-10s %-10s\n",
            "Scenario", "N", "Orig Sens", "Orig Spec", "Orig Acc",
            "Strat Sens", "Strat Spec", "Strat Acc"))
cat(strrep("-", 90), "\n", sep="")

for (sc in 1:6) {
  s <- res[res$sc == sc, ]
  cat(sprintf("%-14s %5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
              sc_names[sc], nrow(s),
              mean(s$orig_sens), mean(s$orig_spec), mean(s$orig_acc),
              mean(s$strat_sens), mean(s$strat_spec), mean(s$strat_acc)))
}

cat("\n")
