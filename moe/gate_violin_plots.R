#!/usr/bin/env Rscript
#' Violin plots: Gate vs CV AUC distributions
#' Generates comparison violin plots for the MoE-K gate paper.
#' Output: moe/results/gate_vs_cv_violin.png

library(ggplot2)
library(data.table)

# Load evaluation data
wk <- if (dir.exists("results")) "." else "moe"
csv_path <- file.path(wk, "results", "gate_evaluation.csv")
if (!file.exists(csv_path)) {
  # Fallback to workspace copy
  alt <- "~/.openclaw/workspace/omen_results/gate_evaluation_2026-06-02.csv"
  csv_path <- path.expand(alt)
}
d <- fread(csv_path)

dt <- melt(d, id.vars = c("family", "n_train"),
           measure.vars = c("gate_auc", "cv_auc"),
           variable.name = "method", value.name = "auc")
dt[, method := fifelse(method == "gate_auc", "Gate", "CV")]
dt[, n_train := factor(n_train, levels = c(200, 300, 500, 1000))]

n_labels <- dt[, .(N = .N / 2), by = n_train]
n_labels[, label := paste0("n = ", n_train, "\n(N = ", N, ")")]

# === Plot 1: Gate vs CV by sample size ===
p1 <- ggplot(dt[!is.na(auc)], aes(x = method, y = auc, fill = method)) +
  geom_violin(trim = TRUE, bw = 0.04, alpha = 0.6) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.5, outlier.size = 0.5) +
  facet_wrap(~ n_train, ncol = 4, labeller = label_both) +
  scale_fill_manual(values = c("Gate" = "#4682B4", "CV" = "#CD5C5C")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0.1, 1.0)) +
  labs(title = "Gate vs Cross-Validation AUC Distribution by Sample Size",
       subtitle = paste0("989 scenarios · Gate beats CV 28.8% vs 19.9%, 51.4% ties"),
       y = "AUC", x = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(wk, "results", "gate_vs_cv_violin.png"),
       p1, width = 12, height = 5, dpi = 150)
cat(sprintf("Saved: results/gate_vs_cv_violin.png\n"))

# === Plot 2: Gate-Oracle gap by family ===
dt_gap <- copy(d)
dt_gap[, gap := gate_auc - oracle_auc]

p2 <- ggplot(dt_gap[!is.na(gap)], aes(x = family, y = gap, fill = family)) +
  geom_violin(trim = TRUE, bw = 0.03, alpha = 0.6) +
  geom_boxplot(width = 0.12, fill = "white", alpha = 0.5, outlier.size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = -0.10, linetype = "dotted", color = "gray60", linewidth = 0.4) +
  annotate("text", x = 0.5, y = -0.095, label = "Δ = −0.10", hjust = 0, vjust = 1,
           size = 3, color = "gray40") +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(-0.4, 0.05), breaks = seq(-0.4, 0, 0.1)) +
  labs(title = "Gate-Oracle AUC Gap by Scenario Family",
       subtitle = paste0("Dashed line: oracle match (0) · Dotted line: −0.10 threshold · ",
                         sum(dt_gap$gap < -0.10, na.rm = TRUE), "/", nrow(dt_gap), " scenarios below"),
       y = "Gate AUC − Oracle AUC", x = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(wk, "results", "gate_gap_by_family.png"),
       p2, width = 10, height = 5.5, dpi = 150)
cat(sprintf("Saved: results/gate_gap_by_family.png\n"))

cat("\nDone. Generated:\n")
cat("  results/gate_vs_cv_violin.png\n")
cat("  results/gate_gap_by_family.png\n")
