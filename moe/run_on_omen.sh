#!/bin/bash
# Run MoE-K simulation on HP Omen (RTX 5090, Windows 11)
#
# 1. Pulls the moe-integration branch from GitHub
# 2. Launches the 1,000-rep simulation pipeline (packages must be pre-installed)
# 3. Monitors progress, checks for errors early
# 4. Pushes results back
#
# Prerequisites (run once in RStudio or RGui):
#   install.packages(c("mvtnorm","glmnet","survival","xgboost","ranger","dplyr","moments"),
#                     repos = "https://cloud.r-project.org")
#
# Run from Git Bash:
#   cd ~/AIProjects/stratified-RMST-boosting && bash moe/run_on_omen.sh

set -e

REPO_DIR="$HOME/AIProjects/stratified-RMST-boosting"
BRANCH="moe-integration"
RSCRIPT="/c/Program Files/R/R-4.6.0/bin/Rscript.exe"

echo "=== MoE-K Simulation on Omen ==="
echo "Time: $(date)"
echo "Branch: $BRANCH"
echo "Repo: $REPO_DIR"
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

# Step 3: Run the simulation
echo ""
echo "=== Starting simulation ==="
echo "Configs: 1,000, Reps per config: 5, Total: 5,000"
echo "Boost iterations (nr): 30"
echo ""

cd "$REPO_DIR"

"$RSCRIPT" -e '
N_CONFIGS <<- 1000
N_REPS <<- 5
source("moe/moe_simulation.R", local=TRUE)
' > moe/simulation_output.log 2>&1 &

PID=$!
echo "Simulation PID: $PID"
echo "Log: moe/simulation_output.log"
echo ""

# Wait a few seconds then check for early failure
sleep 3
if kill -0 $PID 2>/dev/null; then
    echo "R process still running after 3s — looks good ✓"
else
    echo "R process exited early! Log contents:"
    echo "----------------------------------------"
    cat moe/simulation_output.log 2>/dev/null || echo "(log file empty or missing)"
    echo "----------------------------------------"
    echo "Check the log above for R error messages."
    echo "Common issues: missing packages, source() path errors."
    exit 1
fi

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