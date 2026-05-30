#!/bin/bash
# Run MoE-K simulation on HP Omen (RTX 5090, Windows 11)
#
# 1. Pulls the moe-integration branch from GitHub
# 2. Launches the 1,000-rep simulation pipeline (packages must be pre-installed)
# 3. Pushes results back
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

# Step 3: Run the simulation via wrapper (no bash redirects — R handles logging via sink())
echo ""
echo "=== Starting simulation ==="
echo "Configs: 1,000, Reps per config: 5, Total: 5,000"
echo "Boost iterations (nr): 30"
echo ""
echo "Log: moe/simulation_output.log"
echo ""

"$RSCRIPT" moe/run_simulation.R

EXIT_CODE=$?
echo ""
echo "R process exited with code: $EXIT_CODE"

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    cat moe/simulation_output.log 2>/dev/null | tail -50
    echo ""
    echo "=== ERROR: R script failed with exit code $EXIT_CODE ==="
    exit $EXIT_CODE
fi

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