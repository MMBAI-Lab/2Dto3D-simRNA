#!/bin/bash
# Orchestrator: wait for the current example SimRNA to finish, then run the
# three APT-PF1 sequences sequentially (one replica-exchange bundle at a time
# so the server is not saturated). Designed to be launched with `nohup ... &`
# and survive the parent shell exiting.
#
# Log this script's progress to results/APT-PF1/_orchestrator.log.

set -uo pipefail

REPO=/home/pdans/Documents/2Dto3D-simRNA
SIMRNA_DIR="$REPO/SimRNA_64bitIntel_Linux"
export PATH="$SIMRNA_DIR:$PATH"

CONFIG_SRC="$REPO/inputs/example/config.in"
EXAMPLE_PID=2595807
POLL_SECS=60

APT_ROOT_IN="$REPO/inputs/APT-PF1"
APT_ROOT_OUT="$REPO/results/APT-PF1"
mkdir -p "$APT_ROOT_OUT"
LOG="$APT_ROOT_OUT/_orchestrator.log"

log() {
    printf '[%s] %s\n' "$(date -Is)" "$*" >> "$LOG"
}

# 1) Wait for the example run to finish.
log "Orchestrator started. Waiting for example SimRNA (PID=$EXAMPLE_PID)."
while kill -0 "$EXAMPLE_PID" 2>/dev/null; do
    sleep "$POLL_SECS"
done
log "Example run (PID=$EXAMPLE_PID) no longer alive; continuing."

# 2) Define the three APT-PF1 runs (name -> sequence + dot-bracket).
#    Names sanitized for filesystem use (spaces/commas stripped).
NAMES=(
    "DNA_as_RNA_VFold2D"
    "DNA_as_RNA_NUPACK4"
    "NA_as_DNA_NUPACK4"
)
SEQS=(
    "AACACACACAAGGAAGAUCGGACUGAUCCAUAAGGGAAAU"
    "AACACACACAAGGAAGAUCGGACUGAUCCAUAAGGGAAAU"
    "AACACACACAAGGAAGAUCGGACUGAUCCAUAAGGGAAAU"
)
STRS=(
    "..........(((........)))..((((...))))..."
    "...............((((.....))))............"
    "................(((.....)))((.....))...."
)

# 3) Run each in turn.
for i in 0 1 2; do
    name="${NAMES[$i]}"
    seq="${SEQS[$i]}"
    str="${STRS[$i]}"
    indir="$APT_ROOT_IN/$name"
    rundir="$APT_ROOT_OUT/$name"

    log "===== Run $((i+1))/3: $name ====="

    mkdir -p "$indir" "$rundir"

    # Write per-run source inputs under inputs/APT-PF1/<name>/
    printf '%s\n' "$seq" > "$indir/example"
    printf '%s\n' "$str" > "$indir/example.str"
    cp "$CONFIG_SRC" "$indir/config.in"

    # Stage into the run dir (data symlink + copy of inputs).
    ln -sfn "$SIMRNA_DIR/data" "$rundir/data"
    cp "$indir/example" "$indir/example.str" "$indir/config.in" "$rundir/"

    log "  inputs staged in $rundir; launching SimRNA (REMC, 16 replicas)"

    # Launch SimRNA and wait for it to finish before moving on.
    cd "$rundir" || { log "  FATAL: cannot cd into $rundir"; continue; }
    SimRNA -s example -S example.str -c config.in -E 16 -o example >& example.log &
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

log "All APT-PF1 runs attempted. Orchestrator done."
