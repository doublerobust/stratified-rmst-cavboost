#!/usr/bin/env python3
"""Two-phase gate evaluation.

Phase 1: Train RF gate once (sequential).
Phase 2: Split test reps across N cores for real 5-fold CV (parallel).

Usage:
    python3 moe/run_parallel_gate.py              # auto: 10 cores
    python3 moe/run_parallel_gate.py --cores 8    # custom core count
    python3 moe/run_parallel_gate.py --skip-train # skip Phase 1 (reuse saved model)
    python3 moe/run_parallel_gate.py --dry-run    # preview splits only
"""
import argparse
import os
import subprocess
import sys
import time

REPO = subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True
).strip()
TRAIN_SCRIPT = os.path.join(REPO, "moe", "train_gate.R")
EVAL_SCRIPT = os.path.join(REPO, "moe", "evaluate_gate.R")
DEFAULT_CORES = 10


def run_phase1() -> bool:
    """Train the RF gate. Returns True on success."""
    print("=" * 60)
    print("Phase 1: Training RF gate (sequential)")
    print("=" * 60)
    t0 = time.time()
    result = subprocess.run(
        ["Rscript", TRAIN_SCRIPT],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    elapsed = time.time() - t0
    # Print stdout
    for line in result.stdout.splitlines():
        print(f"  {line}")
    if result.stderr:
        for line in result.stderr.splitlines():
            print(f"  [stderr] {line}")
    if result.returncode != 0:
        print(f"\n  ERROR: Phase 1 failed (code {result.returncode}) in {elapsed:.0f}s")
        return False
    print(f"\n  Phase 1 done in {elapsed:.0f}s\n")
    return True


def run_chunk(chunk_idx: int, n_chunks: int, log_dir: str) -> subprocess.Popen:
    """Launch one R process for a chunk. Returns Popen handle."""
    log_file = os.path.join(log_dir, f"chunk_{chunk_idx}.log")
    log_fh = open(log_file, "w")
    proc = subprocess.Popen(
        ["Rscript", EVAL_SCRIPT, str(chunk_idx), str(n_chunks)],
        cwd=REPO,
        stdout=log_fh,
        stderr=subprocess.STDOUT,
    )
    return proc


def merge_results(n_chunks: int) -> str:
    """Merge chunk CSVs into one final CSV. Returns path."""
    import pandas as pd

    chunks = []
    for i in range(n_chunks):
        f = os.path.join(REPO, "moe", "results", f"gate_evaluation_chunk_{i}.csv")
        if os.path.exists(f):
            df = pd.read_csv(f)
            chunks.append(df)
            os.remove(f)
        else:
            print(f"  WARNING: chunk {i} file not found: {f}")

    if not chunks:
        print("  ERROR: no chunk results to merge!")
        return ""

    merged = pd.concat(chunks, ignore_index=True)
    out = os.path.join(REPO, "moe", "results", "gate_evaluation.csv")
    merged.to_csv(out, index=False)
    print(f"\n  Merged {len(merged)} rows into {out}")
    return out


def main():
    parser = argparse.ArgumentParser(description="Parallel gate evaluation")
    parser.add_argument("--cores", type=int, default=DEFAULT_CORES,
                        help=f"Number of parallel cores (default: {DEFAULT_CORES})")
    parser.add_argument("--skip-train", action="store_true",
                        help="Skip Phase 1 (reuse existing saved model)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show chunk splits without running")
    args = parser.parse_args()

    n_chunks = args.cores
    print(f"Gate evaluation: {n_chunks} parallel chunks")
    print(f"   Train: {TRAIN_SCRIPT}")
    print(f"   Eval:  {EVAL_SCRIPT}")
    print(f"   CWD:   {REPO}")
    print()

    # Create log dir
    log_dir = os.path.join(REPO, "moe", "results", "logs")
    os.makedirs(log_dir, exist_ok=True)

    # Dry run
    if args.dry_run:
        import glob
        files = sorted(glob.glob(os.path.join(REPO, "moe", "results", "rep_*.rds")))
        total = len(files)
        chunk_size = (total + n_chunks - 1) // n_chunks
        print(f"   {total} test files across {n_chunks} chunks:\n")
        for i in range(n_chunks):
            start = i * chunk_size
            end = min((i + 1) * chunk_size, total)
            print(f"   Chunk {i+1:2d}: files {start+1:4d}-{end:4d} ({end-start:3d} reps)")
        print()
        return

    # ---- Phase 1: Train RF gate ----
    if not args.skip_train:
        if not run_phase1():
            sys.exit(1)
    else:
        print("Skipping Phase 1 (--skip-train)")

    # ---- Phase 2: Parallel CV ----
    print("=" * 60)
    print("Phase 2: Real 5-fold CV in parallel chunks")
    print("=" * 60)
    t0 = time.time()

    procs = []
    for i in range(n_chunks):
        p = run_chunk(i, n_chunks, log_dir)
        procs.append(p)

    done = [False] * n_chunks
    while not all(done):
        for i, p in enumerate(procs):
            if not done[i]:
                ret = p.poll()
                if ret is not None:
                    done[i] = True
                    status = "OK" if ret == 0 else f"FAIL (code {ret})"
                    elapsed = time.time() - t0
                    print(f"  [{elapsed:6.1f}s] Chunk {i+1}/{n_chunks}: {status}")
        if not all(done):
            time.sleep(5)

    total_time = time.time() - t0
    print(f"\nAll chunks done in {total_time:.0f}s ({total_time/60:.1f}m)")

    # Merge results
    print("\nMerging chunk results...")
    merged = merge_results(n_chunks)
    if merged:
        print(f"Final results: {merged}")
    else:
        print(f"Merge failed -- check logs in {log_dir}/")

    print(f"\nLogs: {log_dir}/chunk_*.log")


if __name__ == "__main__":
    main()
