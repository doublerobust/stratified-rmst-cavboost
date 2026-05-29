# Sequential runner: runs each rep as a separate R process
# Usage: Rscript run_sequential.R [scenario_to_resume_from]
# Each rep saves to results/sc{scenario}_rep{rep}.rds

args <- commandArgs(trailingOnly = TRUE)
start_sc <- if (length(args) >= 1) as.numeric(args[1]) else 1
dir.create("results", showWarnings = FALSE)

sc_names <- c("1"="S1_Linear","2"="S2_Diff","3"="S3_U",
              "4"="S4_Enclave","5"="S5_S","6"="S6_Cross")

for (sc in start_sc:6) {
  sn <- sc_names[[as.character(sc)]]
  done <- length(list.files("results", sprintf("sc%d_rep.*\\.rds", sc)))
  if (done >= 50) { cat(sn, "already done\n"); next }
  
  for (rep in 1:50) {
    out <- sprintf("results/sc%d_rep%d.rds", sc, rep)
    if (file.exists(out)) next  # skip completed reps
    
    code <- system2("Rscript", c("run_one_rep_v2.R", sc, rep, out), stdout = FALSE, stderr = FALSE)
    if (code != 0) {
      cat(sprintf("FAILED: %s rep %d (code %d)\n", sn, rep, code))
      Sys.sleep(5)  # brief pause before retry
    }
    cat(sprintf("\r%s rep %d/50", sn, rep))
  }
  cat(sprintf("\n%s done\n", sn))
}

# Summary
library(data.table)
cat("\n=====================================\n")
cat("  MAIN COMPARISON (50 reps each)\n")
cat("=====================================\n\n")
for (sc in 1:6) {
  sn <- sc_names[[as.character(sc)]]
  files <- list.files("results", sprintf("sc%d_rep.*\\.rds", sc), full.names = TRUE)
  if (length(files) == 0) { cat(sprintf("%s: no results\n\n", sn)); next }
  res <- do.call(rbind, lapply(files, function(f) { d <- readRDS(f); if(is.null(d)) NULL else d }))
  cat(sprintf("--- %s (n=%d) ---\n", sn, nrow(res)))
  cat(sprintf("%-10s %8s\n", "Method", "AUC"))
  for (mn in list(c("orig","Original"), c("strat","Strat_CF"), c("vt","VT"))) {
    v <- mean(res[[paste0(mn[[1]],"_auc")]], na.rm = TRUE)
    cat(sprintf("%-10s %8.4f\n", mn[[2]], v))
  }
  cat("\n")
}
cat("All results saved in results/ directory.\n")
