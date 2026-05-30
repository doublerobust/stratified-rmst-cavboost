#!/usr/bin/env Rscript
# MoE-K quick batch: 50 configs × 1 rep, nr=30
N_CONFIGS <- 50
N_REPS <- 1
NR <- 30

setwd("/home/yue-shentu/workspace/stratified-rmst-cavboost")
source("moe/moe_simulation.R", local = TRUE)
