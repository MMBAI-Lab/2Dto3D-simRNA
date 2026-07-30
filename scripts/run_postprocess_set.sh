#!/bin/bash
# Post-process a whole SimRNA set: for each run, pool the low-T frames, cluster
# by RMSD and by 2D pattern, and build the bilingual branded reports.
#
# Usage: scripts/run_postprocess_set.sh <set_name> [base_name]
#   set_name   Directory name under results/ (e.g. comercialApt).
#   base_name  Prefix of the SimRNA run files (default: example).
#
# Safe to launch while the set is still simulating: each run is picked up as
# soon as scripts/run_simrna_set.sh logs its completion, so finished runs are
# post-processed while later ones are still sampling. Work is sequential — one
# core at a time — so it does not contend with the 16-thread REMC bundle in
# flight.
#
# Per-run work is delegated to scripts/run_cluster_populations.sh, which is
# idempotent: re-running skips whatever already exists.
#
# Launch detached:
#   nohup setsid scripts/run_postprocess_set.sh comercialApt >/dev/null 2>&1 &
#
# Progress goes to results/<set>/_postprocess.log.

set -uo pipefail

SET_NAME="${1:?missing set name (e.g. comercialApt)}"
BASE_NAME="${2:-example}"

REPO="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
OUT_ROOT="$REPO/results/$SET_NAME"
IN_ROOT="$REPO/inputs/$SET_NAME"
ORCH_LOG="$OUT_ROOT/_orchestrator.log"
LOG="$OUT_ROOT/_postprocess.log"
POLL_SECS=120

[ -d "$OUT_ROOT" ] || { echo "no such result set: $OUT_ROOT" >&2; exit 1; }

log() {
    printf '[%s] %s\n' "$(date -Is)" "$*" >> "$LOG"
}

orchestrator_alive() {
    pgrep -f "run_simrna_set.sh $SET_NAME" >/dev/null 2>&1
}

# A run is finished once the orchestrator has logged its exit for that name.
run_finished() {
    grep -qE "  $1 (completed|exited)" "$ORCH_LOG" 2>/dev/null
}

# Same order the orchestrator uses: shortest sequence first, so the earliest
# reports land soonest.
mapfile -t RUNS < <(
    for d in "$IN_ROOT"/*/; do
        [ -f "$d/$BASE_NAME" ] || continue
        printf '%s\t%s\n' "$(tr -d '\n' < "$d/$BASE_NAME" | wc -c)" "$(basename "$d")"
    done | sort -n | cut -f2
)

[ "${#RUNS[@]}" -gt 0 ] || { echo "no staged runs found under $IN_ROOT" >&2; exit 1; }

log "===== Post-processing set '$SET_NAME': ${#RUNS[@]} run(s)"
log "Order: ${RUNS[*]}"

for i in "${!RUNS[@]}"; do
    name="${RUNS[$i]}"
    rundir="$OUT_ROOT/$name"

    log "===== Run $((i+1))/${#RUNS[@]}: $name ====="

    # Wait for this run's simulation to be done.
    waited=0
    while ! run_finished "$name"; do
        if ! orchestrator_alive; then
            log "  SKIP: orchestrator is gone and '$name' has no completion line in $(basename "$ORCH_LOG")"
            continue 2
        fi
        [ "$waited" -eq 0 ] && log "  waiting for the simulation to finish..."
        waited=1
        sleep "$POLL_SECS"
    done
    [ "$waited" -eq 1 ] && log "  simulation finished; starting post-processing"

    if [ ! -s "$rundir/$BASE_NAME.log" ]; then
        log "  SKIP: $rundir/$BASE_NAME.log missing or empty"
        continue
    fi

    start=$(date +%s)
    "$REPO/scripts/run_cluster_populations.sh" "$rundir" "$BASE_NAME" >> "$LOG" 2>&1
    rc=$?
    mins=$(( ($(date +%s) - start) / 60 ))

    if [ "$rc" -eq 0 ]; then
        log "  $name post-processed OK in ${mins} min"
        for lang in EN ES; do
            [ -s "$rundir/README_${lang}.md" ] || log "  WARNING: $rundir/README_${lang}.md was not written"
        done
    else
        log "  $name post-processing FAILED (exit=$rc) after ${mins} min"
    fi
done

log "All '$SET_NAME' runs post-processed. Reports: results/$SET_NAME/*/README_{EN,ES}.md"
