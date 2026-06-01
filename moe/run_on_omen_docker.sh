#!/bin/bash
# run_on_omen_docker.sh — MoE-K simulation via Docker on Omen
#
# Prerequisites:
#   1. Install Docker Desktop for Windows from docker.com
#   2. On agent-server: docker build -t moe-k-sim moe/
#   3. docker push moe-k-sim (or save/load the image on Omen)
#
# Run from Git Bash or PowerShell:
#   bash moe/run_on_omen_docker.sh

set -e

REPO_DIR="$HOME/AIProjects/stratified-RMST-boosting"
RESULTS_DIR="$REPO_DIR/moe/results"
IMAGE="moe-k-sim"      # local build, or "ghcr.io/doublerobust/moe-k-sim" if pushed

echo "=== MoE-K Simulation via Docker ==="
echo "Time: $(date)"
echo "Results: $RESULTS_DIR"
echo ""

# Ensure repo is up to date
if [ -d "$REPO_DIR" ]; then
    cd "$REPO_DIR"
    git stash
    git checkout moe-integration
    git pull origin moe-integration
    echo "Repo updated ✅"
else
    git clone --branch moe-integration \
        https://github.com/doublerobust/stratified-rmst-cavboost.git "$REPO_DIR"
    echo "Repo cloned ✅"
fi

mkdir -p "$RESULTS_DIR"

echo ""
echo "=== Running simulation ($N_WORKERS workers) ==="
echo ""

# Run the gate evaluation (parallel chunks within container)
docker run --rm \
    -v "$REPO_DIR":/app \
    -v "$RESULTS_DIR":/app/moe/results \
    --cpus 11 \
    --memory 16g \
    "$IMAGE" \
    Rscript moe/evaluate_gate.R 0 10

echo ""
echo "=== Running gate importance visualization ==="
echo ""

docker run --rm \
    -v "$REPO_DIR":/app \
    -v "$RESULTS_DIR":/app/moe/results \
    "$IMAGE" \
    Rscript moe/gate_importance_viz.R

echo ""
echo "=== All done ==="
echo "Results saved to: $RESULTS_DIR"
