#!/bin/bash
# build.sh — Build and push MoE-K Docker image
# Run on agent-server (Linux). Pull on Omen (Windows).
#
# Usage:
#   bash moe/build.sh               # build locally
#   bash moe/build.sh --push        # build + push to ghcr.io

set -e

IMAGE_NAME="moe-k-sim"
REGISTRY="ghcr.io/doublerobust"
IMAGE_TAG="${REGISTRY}/${IMAGE_NAME}:latest"

cd "$(dirname "$0")/.."

echo "=== Building Docker image ==="
docker build -t "$IMAGE_NAME" -f moe/Dockerfile .

if [[ "$1" == "--push" ]]; then
    echo "=== Tagging and pushing to $IMAGE_TAG ==="
    docker tag "$IMAGE_NAME" "$IMAGE_TAG"
    docker push "$IMAGE_TAG"
    echo "Pushed ✅"
fi

echo ""
echo "Local image built: $IMAGE_NAME"
echo "To use: docker run --rm -v ... moe-k-sim Rscript moe/evaluate_gate.R"
