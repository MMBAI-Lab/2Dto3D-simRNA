#!/bin/bash
# Run a whole set of SimRNA REMC simulations sequentially.
#
# A "set" is a directory under inputs/ holding one subdirectory per run, each
# staged by scripts/prepare_runs_from_efasta.py with `example`, `example.str`
# and `config.in`. This generalizes the one-off APT-PF1 orchestrator.
#
# Usage: scripts/run_simrna_set.sh <set_name> [run_name ...]
#   set_name   Directory name under inputs/ and results/ (e.g. comercialApt).
#   run_name   Optional explicit run order. Defaults to shortest sequence
#              first, so a staging mistake surfaces in minutes rather than
#              after the longest job.
#
# Runs are sequential on purpose: SimRNA is OpenMP-threaded and a 16-replica
# bundle already saturates most of this box. Two concurrent bundles oversubscribe
# the CPUs and both runs get slower.
#
# Launch detached:
#   nohup scripts/run_simrna_set.sh comercialApt >/dev/null 2>&1 &
#
# Progress goes to results/<set>/_orchestrator.log.

set -uo pipefail

SET_NAME="${1:?missing set name (e.g. comercialApt)}"
shift || true

REPO="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
SIMRNA_DIR="$REPO/SimRNA_64bitIntel_Linux"
export PATH="$SIMRNA_DIR:$PATH"

IN_ROOT="$REPO/inputs/$SET_NAME"
OUT_ROOT="$REPO/results/$SET_NAME"
REPLICAS=16

[ -d "$IN_ROOT" ] || { echo "no such input set: $IN_ROOT" >&2; exit 1; }
mkdir -p "$OUT_ROOT"
LOG="$OUT_ROOT/_orchestrator.log"

log() {
    printf '[%s] %s\n' "$(date -Is)" "$*" >> "$LOG"
}

# Determine run order: explicit args, else shortest sequence first.
if [ "$#" -gt 0 ]; then
    RUNS=("$@")
else
    mapfile -t RUNS < <(
        for d in "$IN_ROOT"/*/; do
            [ -f "$d/example" ] || continue
            printf '%s\t%s\n' "$(tr -d '\n' < "$d/example" | wc -c)" "$(basename "$d")"
        done | sort -n | cut -f2
    )
fi

[ "${#RUNS[@]}" -gt 0 ] || { echo "no staged runs found under $IN_ROOT" >&2; exit 1; }

log "===== Set '$SET_NAME': ${#RUNS[@]} run(s), sequential, $REPLICAS replicas each"
log "Order: ${RUNS[*]}"

for i in "${!RUNS[@]}"; do
    name="${RUNS[$i]}"
    indir="$IN_ROOT/$name"
    rundir="$OUT_ROOT/$name"

    log "===== Run $((i+1))/${#RUNS[@]}: $name ====="

    if [ ! -f "$indir/example" ] || [ ! -f "$indir/example.str" ] || [ ! -f "$indir/config.in" ]; then
        log "  SKIP: $indir is missing example / example.str / config.in"
        continue
    fi
    if [ -s "$rundir/example.log" ]; then
        log "  SKIP: $rundir/example.log already exists (run already done — delete it to redo)"
        continue
    fi

    mkdir -p "$rundir"
    # Every SimRNA binary refuses to start without `data/` in the CWD.
    ln -sfn "$SIMRNA_DIR/data" "$rundir/data"
    cp "$indir/example" "$indir/example.str" "$indir/config.in" "$rundir/"

    nt=$(tr -d '\n' < "$indir/example" | wc -c)
    log "  ${nt} nt · 2D $(cat "$indir/example.str")"
    log "  staged in $rundir; launching SimRNA (REMC, $REPLICAS replicas)"

    cd "$rundir" || { log "  FATAL: cannot cd into $rundir"; continue; }
    SimRNA -s example -S example.str -c config.in -E "$REPLICAS" -o example >& example.log &
    pid=$!
    log "  PID=$pid"

    wait "$pid"
    rc=$?
    if [ "$rc" -eq 0 ]; then
        log "  $name completed cleanly (exit=$rc)"
    else
        log "  $name exited with code $rc (continuing to next run)"
    fi
done

log "All '$SET_NAME' runs attempted. Orchestrator done."
