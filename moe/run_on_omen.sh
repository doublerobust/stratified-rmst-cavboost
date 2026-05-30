#!/bin/bash
# Run MoE-K simulation on HP Omen (RTX 5090, Windows 11)
#
# Launches 11 parallel R processes, each processing a slice of configs.
# No PSOCK cluster needed — plain OS-level multiprocessing.
#
# Prerequisites (run once in RStudio or RGui):
#   install.packages(c("mvtnorm","glmnet","survival","xgboost","ranger","dplyr","moments"),
#                     repos = "https://cloud.r-project.org")
#
# Run from Git Bash:
#   cd ~/AIProjects/stratified-RMST-boosting && bash moe/run_on_omen.sh

# No set -e — we handle errors explicitly

REPO_DIR="$HOME/AIProjects/stratified-RMST-boosting"
BRANCH="moe-integration"
RSCRIPT="/c/Program Files/R/R-4.6.0/bin/Rscript.exe"
N_WORKERS=11

echo "=== MoE-K Simulation on Omen ==="
echo "Time: $(date)"
echo "Branch: $BRANCH"
echo "Repo: $REPO_DIR"
echo "Workers: $N_WORKERS"
echo ""

# Step 1: Clone or pull
if [ -d "$REPO_DIR" ]; then
    cd "$REPO_DIR"
    git stash
    git checkout "$BRANCH"
    git pull origin "$BRANCH"
    echo "Repo updated ✅"
else
    git clone --branch "$BRANCH" https://github.com/doublerobust/stratified-rmst-cavboost.git "$REPO_DIR"
    cd "$REPO_DIR"
    echo "Repo cloned ✅"
fi

# Step 2: Create output dirs
mkdir -p moe/raw moe/results

# Step 3: Launch parallel workers
echo ""
echo "=== Starting simulation ==="
echo "Configs: 1,000, Reps per config: 5, Total: 5,000"
echo "Workers: $N_WORKERS parallel R processes"
echo "Boost iterations (nr): 30"
echo ""

# Calculate config ranges for each worker
# 1000 / 11 = ~91 per worker
CONFIGS_PER_WORKER=$((1000 / N_WORKERS))
PIDS=""

for i in $(seq 1 $N_WORKERS); do
  START=$(((i - 1) * CONFIGS_PER_WORKER + 1))
  if [ $i -eq $N_WORKERS ]; then
    END=1000
  else
    END=$((i * CONFIGS_PER_WORKER))
  fi
  LOG_SUFFIX=$i
  
  echo "  Worker $i: configs $START..$END -> log simulation_output_${LOG_SUFFIX}.log"
  
  "$RSCRIPT" moe/run_simulation.R $START $END --log=$LOG_SUFFIX &
  PIDS="$PIDS $!"
done

echo ""
echo "All $N_WORKERS workers launched. Monitoring..."

# Monitor: check every 30s how many workers are still running
while true; do
  RUNNING=""
  for pid in $PIDS; do
    if kill -0 $pid 2>/dev/null; then
      RUNNING="$RUNNING $pid"
    fi
  done
  
  RUNNING_COUNT=$(echo $RUNNING | wc -w | tr -d ' ')
  if [ "$RUNNING_COUNT" -eq 0 ]; then
    break
  fi
  
  echo "[$(date)] $RUNNING_COUNT / $N_WORKERS workers still running"
  
  # Show latest line from each worker's log
  for i in $(seq 1 $N_WORKERS); do
    LOG="moe/simulation_output_${i}.log"
    if [ -f "$LOG" ]; then
      LAST=$(tail -1 "$LOG" 2>/dev/null)
      if [ -n "$LAST" ]; then
        echo "  Worker $i: $LAST"
      fi
    fi
  done
  
  sleep 30
done

echo ""
echo "=== All workers completed ==="

# Check for errors
HAS_ERROR=0
for i in $(seq 1 $N_WORKERS); do
  LOG="moe/simulation_output_${i}.log"
  if grep -q "ERROR\|Error\|error" "$LOG" 2>/dev/null; then
    echo "⚠️  Errors detected in worker $i log ($LOG)"
    grep "ERROR\|Error" "$LOG" 2>/dev/null | tail -5
    HAS_ERROR=1
  fi
done

# Show summary counts
TOTAL_FILES=$(ls -1 moe/results/rep_*.rds 2>/dev/null | wc -l | tr -d ' ')
echo "Total result files: $TOTAL_FILES / 5000"
echo ""

# Step 4: Run final summary (uses all results regardless of which worker produced them)
echo "=== Generating summary ==="
"$RSCRIPT" -e '
RESULTS_DIR <- "moe/results"
result_files <- list.files(RESULTS_DIR, "rep_.*\\.rds$", full.names = TRUE)
if (length(result_files) > 0) {
  summary_list <- lapply(result_files, function(f) {
    r <- readRDS(f)
    data.frame(
      family = r$config$family, n_train = r$config$n_train,
      auc_K1 = r$aucs[1], auc_K2 = r$aucs[2], auc_K3 = r$aucs[3],
      auc_K4 = r$aucs[4], auc_K5 = r$aucs[5],
      oracle_optimal_K = r$oracle_optimal_K, stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, summary_list)
  cat(sprintf("Reps: %d\n", nrow(summary_df)))
  cat("Optimal K distribution:\n")
  print(table(summary_df$oracle_optimal_K))
  write.csv(summary_df, file.path(RESULTS_DIR, "summary.csv"), row.names = FALSE)
  cat("Summary saved to:", file.path(RESULTS_DIR, "summary.csv"), "\n")
} else {
  cat("No result files found.\n")
}
'

# Step 5: Push results back to GitHub
echo ""
echo "=== Pushing results ==="
git add moe/results/ moe/simulation_output_*.log
git commit -m "MoE-K simulation results ($(date +%Y-%m-%d))"
git push origin "$BRANCH"

echo ""
echo "=== All done ==="
echo "Results pushed to: https://github.com/doublerobust/stratified-rmst-cavboost/tree/$BRANCH/moe/results"