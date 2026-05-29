#!/bin/bash
# Parallel 50-rep × 6-scenario simulation runner.
# Uses background processes with a job limit (default: 6).
# Each rep is a fresh R process → no OOM.
# Output: results_csv/Sc{scenario}_Rep{rep}.csv
#
# Usage: bash run_parallel.sh [n_jobs]
#   n_jobs: max parallel processes (default 6, set to number of CPU cores)

N_JOBS=${1:-6}
SCENARIOS=(1 2 3 4 5 6)
N_REPS=50

mkdir -p results_csv

run_rep() {
    local sc=$1 rep=$2
    local outfile="results_csv/tmp_Sc${sc}_Rep${rep}.rds"
    
    # Skip if already done (CSV equivalent exists)
    if [ -f "results_csv/Sc${sc}_Rep${rep}.csv" ]; then
        return 0
    fi
    
    Rscript run_one_rep_v2.R "$sc" "$rep" "$outfile" 2>/dev/null
    
    if [ -f "$outfile" ]; then
        # Convert RDS to CSV, then delete RDS
        Rscript -e "
            r <- readRDS('$outfile');
            write.csv(r, 'results_csv/Sc${sc}_Rep${rep}.csv', row.names=FALSE)
        " 2>/dev/null
        rm -f "$outfile"
        echo "  ✅ Sc${sc} Rep${rep} done"
    else
        echo "  ❌ Sc${sc} Rep${rep} FAILED"
    fi
}

echo "Running ${N_REPS} reps × ${#SCENARIOS[@]} scenarios (${N_JOBS} parallel jobs)"
echo "Started: $(date)"
echo ""

# Track running jobs
running=0
for sc in "${SCENARIOS[@]}"; do
    for rep in $(seq 1 $N_REPS); do
        # Wait if at job limit
        while [ "$(jobs -r | wc -l)" -ge "$N_JOBS" ]; do
            sleep 2
        done
        
        # Launch in background
        run_rep "$sc" "$rep" &
    done
done

# Wait for remaining jobs
wait

echo ""
echo "Completed: $(date)"

# Summary
echo ""
echo "=== Summary ==="
for sc in "${SCENARIOS[@]}"; do
    sc_names=("" "S1_Linear" "S2_Diff" "S3_U" "S4_Enclave" "S5_S" "S6_Cross")
    sn="${sc_names[$sc]}"
    n_done=$(ls results_csv/Sc${sc}_Rep*.csv 2>/dev/null | wc -l)
    echo "  ${sn}: ${n_done}/${N_REPS} reps"
done

echo ""
echo "All results in results_csv/"
echo ""
echo "To add more reps later, just increase N_REPS in this script and re-run."
echo "Existing CSV files will be skipped."
