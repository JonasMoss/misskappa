#!/usr/bin/env bash
# Start the alpha three-sim overnight run in a detached tmux session.
# Resumable: skips chunks that already have a summary.csv, so it continues from
# the banked Sim A (and survives a crash/reboot — just rerun this script).
#
#   bash start_overnight.sh           # launch
#   tmux attach -t alphasim           # watch it live (Ctrl-b d to detach)
#   tmux kill-session -t alphasim     # stop it
set -euo pipefail

SESSION=alphasim
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/results/3sims"
LOG="$HERE/results/3sims_run.log"

# --- knobs (i7-1355U: 2 P-cores [cpu 0-3] + 8 E-cores [cpu 4-11], 15W) --------
# Measured: 8 E-cores ~= 1.25x one P-core (power-limited, very sublinear). We pin
# to the E-cores so the 2 P-cores stay free and the machine stays usable; the
# throughput cost vs grabbing P-cores is only ~30%.
B=1000          # replications per cell (uniform with the banked Sim A)
CORES=8         # one worker per E-core
PIN=4-11        # taskset cpu list = the 8 E-cores; "" disables pinning
NGRID=100,250,500,1000,2000
TIMEOUT=30      # per-fit wall-clock kill (s); anti-hang + Sim C straggler bound
MEM=2500        # per-fit RSS kill (MB); belt-and-suspenders (Louis mem now bounded)
# -----------------------------------------------------------------------------

if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "tmux session '$SESSION' already exists — attach with: tmux attach -t $SESSION"
  exit 1
fi

mkdir -p "$OUT"
PINCMD=""
if [ -n "$PIN" ]; then PINCMD="taskset -c $PIN"; fi
tmux new-session -d -s "$SESSION" \
  "cd '$HERE' && $PINCMD Rscript run_3sims.R $B '$OUT' $CORES $NGRID $TIMEOUT $MEM > '$LOG' 2>&1"

echo "Launched tmux session '$SESSION': B=$B, cores=$CORES (pinned to cpu ${PIN:-all}), n=$NGRID, timeout=${TIMEOUT}s, mem=${MEM}MB/fit"
echo "  watch:   tmux attach -t $SESSION    (Ctrl-b then d to detach)"
echo "  log:     tail -f $LOG"
echo "  results: $OUT"
