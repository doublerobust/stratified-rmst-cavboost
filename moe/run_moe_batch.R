#!/usr/bin/env Rscript
# MoE-K quick batch: 50 configs x 1 rep, nr=30
N_CONFIGS <- 50
N_REPS <- 1
NR <- 30

# Find repo root (script is at REPO/moe/run_moe_batch.R)
repo_dir <- tryCatch(
  dirname(dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))),
  error = function(e) getwd()
)
setwd(repo_dir)
source("moe/moe_simulation.R", local = TRUE)
