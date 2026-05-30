#!/usr/bin/env Rscript
#
# Wrapper script for moe_simulation.R
# Supports config ranges for multi-process parallelism:
#   Rscript moe/run_simulation.R                     # all 1000 configs
#   Rscript moe/run_simulation.R 1 91                # configs 1..91
#   Rscript moe/run_simulation.R 92 182 --log=2     # configs 92..182, log to ..._2.log
#

args <- commandArgs(TRUE)
CONFIG_START <- 1
CONFIG_END <- NA  # NA means all configs
LOG_SUFFIX <- ""

if (length(args) >= 1) CONFIG_START <- as.integer(args[1])
if (length(args) >= 2) CONFIG_END <- as.integer(args[2])
if (length(args) >= 3 && grepl("^--log=", args[3])) {
  LOG_SUFFIX <- sub("^--log=", "", args[3])
}

log_file <- if (LOG_SUFFIX == "") {
  file.path("moe", "simulation_output.log")
} else {
  file.path("moe", sprintf("simulation_output_%s.log", LOG_SUFFIX))
}

# Redirect all output to log file
sink(log_file, split = TRUE)

cat("MoE-K Simulation Runner\n")
cat("=======================\n")
cat("Log file:", normalizePath(log_file), "\n")
cat("Time:", format(Sys.time()), "\n")
cat(sprintf("Config range: %d .. %s\n", CONFIG_START, if (is.na(CONFIG_END)) "ALL" else as.character(CONFIG_END)))
cat("\n")

N_CONFIGS <- 1000
N_REPS <- 5
N_WORKERS <- 1
PARALLEL <- FALSE

cat(sprintf("N_CONFIGS: %d\n", N_CONFIGS))
cat(sprintf("N_REPS: %d\n", N_REPS))
cat("\n")

# Source the main simulation
source("moe/moe_simulation.R")

# Close sink
sink()
cat(sprintf("Simulation complete. Log saved to %s\n", log_file))