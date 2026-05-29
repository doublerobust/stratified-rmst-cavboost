#!/bin/bash
# Fully detached simulation runner.
# Run: setsid bash run_sim_detached.sh &

N_JOBS=6
N_REPS=50
mkdir -p results_csv

exec > /tmp/sim_final.log 2>&1
echo "=== Stratified RMST Simulation ==="
echo "Started: $(date)"
echo "PID: $$"
echo ""

for sc in 1 2 3 4 5 6; do
  for rep in $(seq 1 $N_REPS); do
    outcsv="results_csv/Sc${sc}_Rep${rep}.csv"
    [ -f "$outcsv" ] && continue
    
    while [ "$(jobs -r | wc -l)" -ge "$N_JOBS" ]; do
      sleep 5
    done
    
    (
      outrds="results_csv/tmp_Sc${sc}_Rep${rep}.rds"
      Rscript run_one_rep_v2.R "$sc" "$rep" "$outrds" 2>/dev/null
      if [ -f "$outrds" ]; then
        Rscript -e "r <- readRDS('$outrds'); write.csv(r, '$outcsv', row.names=FALSE)" 2>/dev/null
        rm -f "$outrds"
        echo "$(date +%H:%M) Sc${sc} Rep${rep} done"
      else
        echo "$(date +%H:%M) Sc${sc} Rep${rep} FAILED"
      fi
    ) &
  done
  echo "=== Waiting for Sc${sc} ==="
  wait
  echo "=== Sc${sc} complete ==="
done

echo "=== ALL DONE ==="
echo "Finished: $(date)"
