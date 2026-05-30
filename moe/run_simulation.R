#!/usr/bin/env Rscript
#
# Wrapper script for moe_simulation.R
# Handles logging via sink() to avoid bash redirect issues with R 4.6.0 on Windows.
#

log_file <- file.path("moe", "simulation_output.log")

# Redirect all output to log file
sink(log_file, split = TRUE)

cat("MoE-K Simulation Runner\n")
cat("=======================\n")
cat("Log file:", normalizePath(log_file), "\n")
cat("Time:", format(Sys.time()), "\n")
cat("\n")

N_CONFIGS <- 1000
N_REPS <- 5

cat(sprintf("N_CONFIGS: %d\n", N_CONFIGS))
cat(sprintf("N_REPS: %d\n", N_REPS))
cat("\n")

# Source the main simulation
source("moe/moe_simulation.R")

# Close sink
sink()
cat(sprintf("Simulation complete. Log saved to %s\n", log_file))