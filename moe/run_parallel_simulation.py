#!/usr/bin/env python3
"""Parallel MoE simulation: splits configs across N cores.

Usage:
    python3 moe/run_parallel_simulation.py                    # default: N_CONFIGS=500, 12 cores
    python3 moe/run_parallel_simulation.py --configs 200 --cores 8
    python3 moe/run_parallel_simulation.py --configs 500 --reps 10 --cores 12
    python3 moe/run_parallel_simulation.py --dry-run          # preview splits only

After all chunks finish, run merge + extract + evaluate:
    python3 moe/run_parallel_simulation.py --merge-only
    Rscript moe/extract_gate_data.R
    Rscript moe/evaluate_gate.R
"""
import argparse
import os
import subprocess
import sys
import time

REPO = subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True
).strip()
REPO = REPO.replace("\\", "/")
SCRIPT = os.path.join(REPO, "moe", "moe_simulation.R").replace("\\", "/")
LOG_DIR = os.path.join(REPO, "moe", "results", "logs").replace("\\", "/")

DEFAULT_CONFIGS = 500
DEFAULT_REPS = 10
DEFAULT_CORES = 12


def parse_args():
    parser = argparse.ArgumentParser(description="Parallel MoE simulation")
    parser.add_argument("--configs", type=int, default=DEFAULT_CONFIGS,
                        help=f"Total configs (default: {DEFAULT_CONFIGS})")
    parser.add_argument("--reps", type=int, default=DEFAULT_REPS,
                        help=f"Reps per config (default: {DEFAULT_REPS})")
    parser.add_argument("--cores", type=int, default=DEFAULT_CORES,
                        help=f"Parallel cores (default: {DEFAULT_CORES})")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show splits without running")
    parser.add_argument("--clean", action="store_true",
                        help="Delete old RDS/raw files before starting")
    parser.add_argument("--merge-only", action="store_true",
                        help="Skip simulation, just list result files")
    return parser.parse_args()


def run_chunk(chunk_idx: int, start: int, end: int, configs: int,
              reps: int, cores: int, log_dir: str) -> subprocess.Popen:
    """Launch one R process for a config range."""
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"sim_chunk_{chunk_idx}.log")
    log_fh = open(log_file, "w")

    r_code = (
        f"N_CONFIGS <- {configs}; "
        f"N_REPS <- {reps}; "
        f"PARALLEL <- FALSE; "
        f"CONFIG_START <- {start}; "
        f"CONFIG_END <- {end}; "
        f"source('{SCRIPT}')"
    )

    log_fh.write(f"=== Chunk {chunk_idx}: configs {start}-{end} ===\n")
    log_fh.flush()

    proc = subprocess.Popen(
        ["Rscript", "-e", r_code],
        cwd=REPO,
        stdout=log_fh,
        stderr=subprocess.STDOUT,
    )
    return proc


def main():
    args = parse_args()

    if args.merge_only:
        print(f"Result files in {os.path.join(REPO, 'moe', 'results')}:")
        os.system(f"ls {os.path.join(REPO, 'moe', 'results', 'rep_*.rds')} 2>/dev/null | wc -l")
        return

    n_chunks = min(args.cores, args.configs)
    chunk_size = args.configs // n_chunks
    remainder = args.configs % n_chunks

    # Build chunk ranges (1-indexed, inclusive)
    chunks = []
    ptr = 1
    for i in range(n_chunks):
        extra = 1 if i < remainder else 0
        csize = chunk_size + extra
        chunks.append((ptr, ptr + csize - 1))
        ptr += csize

    total = args.configs * args.reps
    print(f"MoE Simulation: {args.configs} configs x {args.reps} reps = {total} total")
    print(f"   Cores: {n_chunks}  Chunk size: ~{chunk_size} configs each")
    print(f"   CWD:   {REPO}")
    print()

    # Preview
    print(f"{'Chunk':>6}  {'Configs':>10}  {'Est. time':>10}")
    print("-" * 30)
    for i, (s, e) in enumerate(chunks):
        print(f"{i+1:>6}  {s:>5}-{e:<5}  ~{((e-s+1)*args.reps*8)//60:>3}m")
    print()

    if args.dry_run:
        return

    # Optional clean
    if args.clean:
        print("Cleaning old RDS and raw files...")
        for d in ["moe/results", "moe/raw"]:
            for f in os.listdir(os.path.join(REPO, d)):
                if f.endswith(".rds"):
                    os.remove(os.path.join(REPO, d, f))
        print("   Done.\n")

    # Launch all chunks
    t0 = time.time()
    procs = []
    for i, (s, e) in enumerate(chunks):
        p = run_chunk(i, s, e, args.configs, args.reps, n_chunks, LOG_DIR)
        procs.append(p)
        print(f"  Launched chunk {i+1}/{n_chunks} (PID {p.pid})")

    # Monitor progress
    done = [False] * n_chunks
    while not all(done):
        for i, p in enumerate(procs):
            if not done[i]:
                ret = p.poll()
                if ret is not None:
                    done[i] = True
                    elapsed = time.time() - t0
                    status = "OK" if ret == 0 else f"FAIL (code {ret})"
                    print(f"  [{elapsed:6.1f}s] Chunk {i+1}/{n_chunks}: {status}")
        if not all(done):
            time.sleep(10)

    total_time = time.time() - t0
    print(f"\nAll chunks done in {total_time:.0f}s ({total_time/60:.1f}m)")

    # Count results
    res_dir = os.path.join(REPO, "moe", "results")
    n_files = len([f for f in os.listdir(res_dir) if f.startswith("rep_") and f.endswith(".rds")])
    print(f"Result files: {n_files} / {total}")
    print(f"Logs: {LOG_DIR}/sim_chunk_*.log")
    print(f"\nNext steps:")
    print(f"  Rscript moe/extract_gate_data.R")
    print(f"  Rscript moe/evaluate_gate.R")


if __name__ == "__main__":
    main()