#!/usr/bin/env python3
"""Parallel simulation runner — 6 workers, CSV output, skip existing."""
import glob, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

REPO = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(REPO, "results_csv")
os.makedirs(CSV_DIR, exist_ok=True)

N_JOBS = int(sys.argv[1]) if len(sys.argv) > 1 else 6
N_REPS = 50
SCENARIOS = [(s, r) for s in range(1, 7) for r in range(1, N_REPS + 1)]


def run_rep(sc, rep):
    out_csv = os.path.join(CSV_DIR, f"Sc{sc}_Rep{rep}.csv")
    if os.path.exists(out_csv):
        return f"Sc{sc}_Rep{rep} — already exists, skipped"

    out_rds = os.path.join(CSV_DIR, f"tmp_Sc{sc}_Rep{rep}.rds")
    result = subprocess.run(
        ["Rscript", "run_one_rep_v2.R", str(sc), str(rep), out_rds],
        capture_output=True, text=True, cwd=REPO
    )
    if result.returncode != 0:
        return f"Sc{sc}_Rep{rep} FAILED: {result.stderr[:200]}"

    if os.path.exists(out_rds):
        subprocess.run(["Rscript", "-e",
            f"r <- readRDS('{out_rds}'); write.csv(r, '{out_csv}', row.names=FALSE)"],
            capture_output=True, cwd=REPO)
        os.remove(out_rds)

    return f"Sc{sc}_Rep{rep} done"


def main():
    start = time.time()
    done = len(glob.glob(os.path.join(CSV_DIR, "Sc*_Rep*.csv")))
    total = len(SCENARIOS)
    remaining = [(s, r) for s, r in SCENARIOS
                 if not os.path.exists(os.path.join(CSV_DIR, f"Sc{s}_Rep{r}.csv"))]
    print(f"📊 {done}/{total} done, {len(remaining)} remaining ({N_JOBS} workers)")
    print(f"   Started: {time.strftime('%H:%M')}")

    with ThreadPoolExecutor(max_workers=N_JOBS) as pool:
        futures = {pool.submit(run_rep, s, r): (s, r) for s, r in remaining}
        n_done = done
        for f in as_completed(futures):
            n_done += 1
            elapsed = time.time() - start
            rate = n_done / elapsed if elapsed > 0 else 0
            eta_sec = (total - n_done) / rate if rate > 0 else 0
            print(f"  [{n_done}/{total}] {f.result()}  "
                  f"({eta_sec/60:.0f}m remaining)")

    elapsed = time.time() - start
    print(f"\n✅ Complete in {elapsed/60:.1f} minutes")
    print(f"   Results in: {CSV_DIR}")


if __name__ == "__main__":
    main()
