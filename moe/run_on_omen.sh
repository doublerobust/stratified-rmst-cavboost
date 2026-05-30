#!/bin/bash
# Run MoE-K simulation on HP Omen (RTX 5090, Windows 11)
#
# 1. Pulls the moe-integration branch from GitHub
# 2. Runs the full 1,000-rep simulation pipeline
# 3. Pushes results back
#
# Run from Windows Terminal (as admin, or in R session with system permissions):
#   wsl bash moe/run_on_omen.sh
#
# Or from RStudio directly:
#   source("moe/run_on_omen.R")

set -e

REPO_DIR="$HOME/stratified-rmst-cavboost"  # adjust to your WSL path
BRANCH="moe-integration"

echo "=== MoE-K Simulation on Omen ==="
echo "Time: $(date)"
echo "Branch: $BRANCH"
echo "Repo: $REPO_DIR"
echo ""

# Step 0: Clone or pull
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

# Step 1: Create output dirs
mkdir -p moe/raw moe/results

# Step 2: Install any missing R packages
echo ""
echo "Checking R packages..."
Rscript -e '
pkgs_ok <- sapply(c("mvtnorm","glmnet","survival","xgboost","ranger","parallel","moments"), 
                    function(p) requireNamespace(p, quietly=TRUE))
missing <- names(pkgs_ok)[!pkgs_ok]
if (length(missing) > 0) {
  cat("Installing missing packages:", missing, "\n")
  install.packages(missing, repos="https://cloud.r-project.org")
} else {
  cat("All packages found ✅\n")
}
'

# Step 3: Run the simulation
echo ""
echo "=== Starting simulation ==="
echo "Configs: 1,000, Reps per config: 5, Total: 5,000"
echo "Boost iterations (nr): 30"
echo ""

cd "$REPO_DIR"

# Run with explicit nr=30 (edit moe_simulation.R to reduce from 50 to 30)
RSCRIPT="/usr/bin/Rscript"  # or wherever R is on your WSL

$RSCRIPT -e '
N_CONFIGS <<- 1000
N_REPS <<- 5
source("moe/moe_simulation.R", local=TRUE)
' > moe/simulation_output.log 2>&1 &

PID=$!
echo "Simulation PID: $PID"
echo "Log: moe/simulation_output.log"
echo ""

# Monitor progress
while kill -0 $PID 2>/dev/null; do
    tail -3 moe/simulation_output.log 2>/dev/null
    echo "--- $(date) ---"
    sleep 60
done

echo ""
echo "=== Simulation complete ==="
tail -20 moe/simulation_output.log

# Step 4: Push results back to GitHub
echo ""
echo "=== Pushing results ==="
git add moe/results/ moe/simulation_output.log
git commit -m "MoE-K simulation results ($(date +%Y-%m-%d))"
git push origin "$BRANCH"

echo ""
echo "=== All done ==="
echo "Results pushed to: https://github.com/doublerobust/stratified-rmst-cavboost/tree/$BRANCH/moe/results"
