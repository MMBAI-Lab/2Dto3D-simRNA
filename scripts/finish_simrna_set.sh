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

# 2) Optional: centroid PDB/PNG per run x cutoff x cluster. Needs the scientific
#    stack; the set report treats these as optional, so a missing numpy must not
#    take the reports down with it.
if python3 -c 'import numpy, matplotlib, Bio' >/dev/null 2>&1; then
    log "--- extracting global centroids"
    python3 "$REPO/scripts/extract_global_centroids.py" --set "$SET_NAME" >> "$LOG" 2>&1
    log "--- centroids returned $?"
else
    log "--- SKIP centroids: numpy/matplotlib/biopython not importable by $(command -v python3)."
    log "    Reports are unaffected (centroid PNGs are optional). To add them later:"
    log "    pip install --user numpy matplotlib biopython && python3 scripts/extract_global_centroids.py --set $SET_NAME"
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
