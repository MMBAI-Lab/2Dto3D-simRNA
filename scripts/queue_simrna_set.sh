#!/bin/bash
# Queue a whole SimRNA set to start once the box is free, then take it all the
# way to the reports.
#
# Usage: scripts/queue_simrna_set.sh <set_name> [--after <other_set>] [run_name ...]
#   set_name     Set to run, staged under inputs/<set_name>/.
#   --after      Wait until this set's orchestrator has finished before starting.
#   run_name...  Optional explicit run order, forwarded to run_simrna_set.sh.
#
# Why queue instead of just launching: a 16-replica REMC bundle already saturates
# this box, so two concurrent sets oversubscribe the CPUs and both get slower
# (same reason run_simrna_set.sh runs its own runs sequentially). The wait is on
# both the named orchestrator *and* any live SimRNA process, so the gap between
# two runs of the set ahead is not mistaken for the box being free.
#
# Launch detached:
#   nohup setsid scripts/queue_simrna_set.sh ApF9983R --after ApF20053R >/dev/null 2>&1 &
#
# Progress goes to results/<set>/_queue.log.

set -uo pipefail

SET_NAME="${1:?missing set name (e.g. ApF9983R)}"
shift

AFTER=""
if [ "${1:-}" = "--after" ]; then
    AFTER="${2:?--after needs a set name}"
    shift 2
fi
RUNS=("$@")

REPO="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
IN_ROOT="$REPO/inputs/$SET_NAME"
OUT_ROOT="$REPO/results/$SET_NAME"
POLL_SECS=120

[ -d "$IN_ROOT" ] || { echo "no such input set: $IN_ROOT" >&2; exit 1; }
[ -n "$AFTER" ] && [ "$AFTER" = "$SET_NAME" ] && {
    echo "--after $AFTER is this same set; that would never start" >&2; exit 1; }

mkdir -p "$OUT_ROOT"
LOG="$OUT_ROOT/_queue.log"

log() {
    printf '[%s] %s\n' "$(date -Is)" "$*" >> "$LOG"
}

box_busy() {
    [ -n "$AFTER" ] && pgrep -f "run_simrna_set\.sh $AFTER" >/dev/null 2>&1 && return 0
    pgrep -x SimRNA >/dev/null 2>&1 && return 0
    return 1
}

log "===== queued set '$SET_NAME'${AFTER:+, waiting for '$AFTER' to finish}"

waited=0
while box_busy; do
    [ "$waited" -eq 0 ] && log "  box busy; polling every ${POLL_SECS}s"
    waited=1
    sleep "$POLL_SECS"
done
[ "$waited" -eq 1 ] && log "  box free; starting"

# 1) Simulate. Backgrounded so the post-processor can overlap with it.
log "--- launching run_simrna_set.sh $SET_NAME ${RUNS[*]:-}"
"$REPO/scripts/run_simrna_set.sh" "$SET_NAME" ${RUNS[@]+"${RUNS[@]}"} &
sim_pid=$!

# run_postprocess_set.sh gives up if the orchestrator is not yet visible, so let
# it register before starting the chain that watches it.
for _ in $(seq 1 60); do
    pgrep -f "run_simrna_set\.sh $SET_NAME" >/dev/null 2>&1 && break
    sleep 1
done

# 2) Post-process + reports, overlapping the simulation run by run.
log "--- launching finish_simrna_set.sh $SET_NAME"
"$REPO/scripts/finish_simrna_set.sh" "$SET_NAME" &
fin_pid=$!

wait "$sim_pid"; log "  simulation orchestrator exited ($?)"
wait "$fin_pid"; log "  reports chain exited ($?)"

log "===== set '$SET_NAME' done end to end"
