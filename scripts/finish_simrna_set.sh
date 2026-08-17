#!/bin/bash
# Finish a SimRNA set: post-process every run, then build the set-level
# bilingual branded report and its DOCX.
#
# Usage: scripts/finish_simrna_set.sh <set_name> [ref_cutoff] [base_name]
#   set_name    Directory name under results/ (e.g. ApF20053R).
#   ref_cutoff  RMSD cutoff (Å) reported in detail by the set report. Default 12,
#               the value used for the 60-nt comercialApt runs; 8 Å over-splits a
#               long chain and 16 Å collapses it.
#   base_name   Prefix of the SimRNA run files (default: example).
#
# Safe to launch while the set is still simulating: run_postprocess_set.sh picks
# up each run as soon as the orchestrator logs its completion, so the last run's
# clustering starts the moment its simulation ends. Work is single-core and does
# not contend with the 16-thread REMC bundle in flight.
#
# Launch detached:
#   nohup setsid scripts/finish_simrna_set.sh ApF20053R >/dev/null 2>&1 &
#
# Progress goes to results/<set>/_postprocess.log (per-run) and
# results/<set>/_reports.log (set-level).

set -uo pipefail

SET_NAME="${1:?missing set name (e.g. ApF20053R)}"
REF_CUTOFF="${2:-12}"
BASE_NAME="${3:-example}"

REPO="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
OUT_ROOT="$REPO/results/$SET_NAME"

[ -d "$OUT_ROOT" ] || { echo "no such result set: $OUT_ROOT" >&2; exit 1; }
LOG="$OUT_ROOT/_reports.log"

log() {
    printf '[%s] %s\n' "$(date -Is)" "$*" >> "$LOG"
}

log "===== finishing set '$SET_NAME' (ref cutoff ${REF_CUTOFF} A)"

# 1) Per-run: pool low-T frames, cluster by RMSD and by 2D, write README_{EN,ES}.
log "--- post-processing runs (see _postprocess.log)"
"$REPO/scripts/run_postprocess_set.sh" "$SET_NAME" "$BASE_NAME"
log "--- post-processing returned $?"

# 2) Centroid PDBs per run x cutoff x cluster, then their renders.
#
#    These are two separate concerns and used to be gated as one, which cost the
#    whole ApF series its centroids: the guard tested numpy/matplotlib/biopython
#    and skipped everything when they were missing. But back-mapping needs only
#    SimRNA_trafl2pdbs -- no Python stack at all -- and the centroid PDBs are the
#    input to the next stage (substate selection for MD), not decoration. So the
#    extraction always runs, and only the drawing is conditional.
log "--- extracting global centroids (PDBs; renders handled separately below)"
python3 "$REPO/scripts/extract_global_centroids.py" --set "$SET_NAME" --no-render >> "$LOG" 2>&1
log "--- centroids returned $?"

#    Rendering, in order of preference. ChimeraX needs nothing from Python and
#    aligns per aptamer with a shared camera, so it is the better output as well
#    as the more available one; the matplotlib tracer is the fallback.
CHIMERAX_BIN="${CHIMERAX_BIN:-/home/pdans/chimerax/bin/chimerax-headless}"
if [ -x "$CHIMERAX_BIN" ]; then
    log "--- rendering centroids with ChimeraX"
    python3 "$REPO/scripts/render_centroids_chimerax.py" --set "$SET_NAME" \
        --chimerax "$CHIMERAX_BIN" >> "$LOG" 2>&1
    log "--- ChimeraX render returned $?"
elif python3 -c 'import numpy, matplotlib, Bio' >/dev/null 2>&1; then
    log "--- ChimeraX not at $CHIMERAX_BIN; rendering with the matplotlib tracer"
    python3 "$REPO/scripts/extract_global_centroids.py" --set "$SET_NAME" >> "$LOG" 2>&1
    log "--- matplotlib render returned $?"
else
    log "--- SKIP renders: no ChimeraX at $CHIMERAX_BIN, and numpy/matplotlib/biopython"
    log "    are not importable by $(command -v python3). The centroid PDBs are still"
    log "    in place and the reports are unaffected -- only the gallery is empty."
    log "    To add it later:  python3 scripts/render_centroids_chimerax.py --set $SET_NAME"
fi

# 3) Set-level cross-run report, bilingual pair + DOCX.
log "--- build_set_report.py ($SET_NAME)"
python3 "$REPO/scripts/build_set_report.py" "$SET_NAME" \
    --default-ref-cutoff "$REF_CUTOFF" >> "$LOG" 2>&1
rc=$?
log "--- build_set_report returned $rc"

if [ "$rc" -eq 0 ]; then
    log "--- build_report.py (DOCX for results/$SET_NAME)"
    python3 "$REPO/scripts/build_report.py" "results/$SET_NAME" >> "$LOG" 2>&1
    log "--- build_report returned $?"
fi

for lang in EN ES; do
    [ -s "$OUT_ROOT/REPORT_${lang}.md" ] || log "WARNING: REPORT_${lang}.md was not written"
done

log "===== set '$SET_NAME' finished. Reports: results/$SET_NAME/{*/README,REPORT}_{EN,ES}.md"
