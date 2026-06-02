#!/usr/bin/env python3
"""Parallel feature extraction: splits RDS files across N cores.

Usage:
    python3 moe/run_parallel_extraction.py                  # auto: 10 cores
    python3 moe/run_parallel_extraction.py --cores 8        # custom
    python3 moe/run_parallel_extraction.py --dry-run        # preview splits

Each chunk writes a partial CSV to moe/results/chunk_N.csv.
The main process merges them into gate_training_data.csv.
"""
import argparse
import os
import subprocess
import sys
import time
import glob

REPO = subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True
).strip()
REPO = REPO.replace("\\", "/")
MOE = os.path.join(REPO, "moe").replace("\\", "/")
RESULTS = os.path.join(MOE, "results").replace("\\", "/")
EXTRACT_SCRIPT = os.path.join(MOE, "extract_gate_data.R").replace("\\", "/")
DEFAULT_CORES = 10


def find_rds_files() -> list:
    """Get sorted list of results RDS files."""
    pattern = os.path.join(RESULTS, "rep_*.rds")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"No RDS files found in {RESULTS}")
        sys.exit(1)
    return files


def split_files(files: list, n_chunks: int) -> list:
    """Split file list into n_chunks roughly equal parts."""
    chunks = [[] for _ in range(n_chunks)]
    for i, f in enumerate(files):
        chunks[i % n_chunks].append(f)
    return chunks


def write_chunk_file(chunk_idx: int, files: list, chunk_dir: str) -> str:
    """Write a text file listing the RDS files for this chunk."""
    os.makedirs(chunk_dir, exist_ok=True)
    path = os.path.join(chunk_dir, f"chunk_{chunk_idx}.txt")
    with open(path, "w") as f:
        for fp in files:
            f.write(fp + "\n")
    return path


def run_chunk(chunk_idx: int, chunk_file: str, log_dir: str) -> subprocess.Popen:
    """Launch one R process for a chunk."""
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"extract_chunk_{chunk_idx}.log")
    log_fh = open(log_file, "w")
    
    # Source the extract script directly — it auto-detects chunk mode from commandArgs
    r_script = f'source("{EXTRACT_SCRIPT}")'

    proc = subprocess.Popen(
        ["Rscript", "-e", r_script, chunk_file, str(chunk_idx)],
        stdout=log_fh,
        stderr=subprocess.STDOUT,
    )
    return proc


def main():
    parser = argparse.ArgumentParser(description="Parallel feature extraction")
    parser.add_argument("--cores", type=int, default=DEFAULT_CORES,
                        help="Number of parallel R processes (default: %d)" % DEFAULT_CORES)
    parser.add_argument("--dry-run", action="store_true",
                        help="Preview splits without running")
    parser.add_argument("--output", default="gate_training_data.csv",
                        help="Output CSV name (default: gate_training_data.csv)")
    args = parser.parse_args()

    files = find_rds_files()
    n_files = len(files)
    n_chunks = min(args.cores, n_files)

    print(f"Found {n_files} RDS files")
    print(f"Splitting across {n_chunks} cores ({n_files // n_chunks}–{n_files // n_chunks + 1} files each)")
    print()

    chunks = split_files(files, n_chunks)

    if args.dry_run:
        for i, chunk in enumerate(chunks):
            print(f"  Chunk {i}: {len(chunk)} files ({chunk[0]} .. {chunk[-1]})")
        return

    chunk_dir = os.path.join(RESULTS, ".extract_chunks")
    log_dir = os.path.join(RESULTS, ".extract_logs")

    # Write chunk files
    chunk_files = []
    for i, chunk in enumerate(chunks):
        cf = write_chunk_file(i, chunk, chunk_dir)
        chunk_files.append(cf)

    # Launch all chunks
    procs = []
    for i in range(n_chunks):
        proc = run_chunk(i, chunk_files[i], log_dir)
        procs.append(proc)
        print(f"  Chunk {i} started (PID {proc.pid}, {len(chunks[i])} files)")

    # Wait for completion
    start = time.time()
    while True:
        done = sum(1 for p in procs if p.poll() is not None)
        elapsed = time.time() - start
        sys.stdout.write(f"\r  Progress: {done}/{n_chunks} chunks complete ({elapsed:.0f}s)")
        sys.stdout.flush()
        if done == n_chunks:
            break
        time.sleep(5)
    print()

    # Check for errors
    errors = [i for i, p in enumerate(procs) if p.returncode != 0]
    if errors:
        print(f"\n❌ {len(errors)} chunks failed: {errors}")
        for i in errors:
            log_path = os.path.join(log_dir, f"extract_chunk_{i}.log")
            print(f"  Check {log_path}")
        sys.exit(1)

    # Merge chunks
    print(f"\nMerging {n_chunks} chunk CSVs...")
    header_written = False
    output_path = os.path.join(RESULTS, args.output)
    with open(output_path, "w") as out:
        for i in range(n_chunks):
            chunk_csv = os.path.join(RESULTS, f"chunk_{i}.csv")
            if not os.path.exists(chunk_csv):
                print(f"  ⚠ Chunk {i} CSV not found, skipping")
                continue
            with open(chunk_csv) as f:
                lines = f.readlines()
                if not header_written:
                    out.write(lines[0])  # header
                    out.writelines(lines[1:])  # data
                    header_written = True
                else:
                    out.writelines(lines[1:])  # skip header
            n_lines = len(lines) - 1
            print(f"  Chunk {i}: {n_lines} rows")
            os.remove(chunk_csv)

    total = sum(1 for _ in open(output_path)) - 1  # minus header
    elapsed = time.time() - start
    print(f"\n✅ Merged: {output_path} ({total} rows, {elapsed:.0f}s)")

    # Cleanup
    for f in chunk_files:
        os.remove(f)
    for d in [chunk_dir, log_dir]:
        if os.path.exists(d):
            os.rmdir(d)

    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
