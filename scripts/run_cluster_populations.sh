#!/bin/bash
# Build per-run conformational substate populations for a completed SimRNA REMC
# run: pool low-T frames → RMSD clustering at four cutoffs → SS-pattern clustering
# → branded README_EN.md + README_ES.md with two populations tables.
#
# Usage: scripts/run_cluster_populations.sh <run_dir> [base_name]
#   run_dir    Directory containing <name>.log, <name>_{01..16}.trafl,
#              <name>_01-000001.pdb, <name>.str, and a `data/` symlink (or this
#              script will create one).
#   base_name  Prefix of the SimRNA run files (default: example).
#
# Re-running is idempotent: skips work whose outputs already exist.

set -euo pipefail

RUN_DIR="$(readlink -f "${1:?missing run_dir}")"
NAME="${2:-example}"

REPO="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
SIMRNA_DIR="$REPO/SimRNA_64bitIntel_Linux"
PKG="$REPO/scripts/simRNA-scripts/scripts"
export PATH="$SIMRNA_DIR:$PATH"

CUTOFFS=(8 12 16)
LOW_T_MIN=2
LOW_T_MAX=4
CLUST_FRACTION=1.0

cd "$RUN_DIR"
[ -e data ] || ln -sfn "$SIMRNA_DIR/data" data

echo "=== $(date -Is) $RUN_DIR (name=$NAME)"

# 1) Pool low-T frames (REMC levels 2-4 per protocol)
if [ ! -s low_temp.trafl ]; then
    echo "--- extract_low_temp_frames: levels ${LOW_T_MIN}-${LOW_T_MAX}"
    python3 "$PKG/extract_low_temp_frames.py" "$NAME.log" \
        --base-name "$NAME" \
        --min-temp "$LOW_T_MIN" --max-temp "$LOW_T_MAX" \
        --output low_temp.trafl
else
    echo "--- low_temp.trafl already exists, skipping pooling"
fi

# 2) RMSD clustering at 4 cutoffs
RMSD_DIR=rmsd_clusters
mkdir -p "$RMSD_DIR"
if [ ! -s "$RMSD_DIR/clustering.log" ]; then
    echo "--- RMSD clustering: cutoffs ${CUTOFFS[*]} Å, fraction $CLUST_FRACTION"
    (
        cd "$RMSD_DIR"
        ln -sfn "$SIMRNA_DIR/data" data
        ln -sfn "../low_temp.trafl" low_temp.trafl
        clustering low_temp.trafl "$CLUST_FRACTION" "${CUTOFFS[@]}" > clustering.log 2>&1
    )
else
    echo "--- rmsd_clusters/clustering.log already present, skipping clustering"
fi

# 3) SS-pattern clustering: rebuild per-frame .ss_detected via SimRNA_trafl2pdbs
SS_DIR=ss_clusters
mkdir -p "$SS_DIR"
if [ ! -s "$SS_DIR/ss_clusters.tsv" ]; then
    echo "--- SS clustering: generating per-frame .ss_detected (this is the slow step)"
    (
        cd "$SS_DIR"
        ln -sfn "$SIMRNA_DIR/data" data
        ln -sfn "../low_temp.trafl" low_temp.trafl
        cp -f "../${NAME}_01-000001.pdb" ref.pdb
        if ! ls ./*.ss_detected >/dev/null 2>&1; then
            SimRNA_trafl2pdbs ref.pdb low_temp.trafl : > trafl2pdbs.log 2>&1 || true
        fi
        python3 "$REPO/scripts/ss_cluster.py" . > ss_clusters.tsv
        # Keep only the .ss_detected files and the TSV; drop the PDB flood.
        find . -maxdepth 1 -name '*.pdb' ! -name ref.pdb -delete
    )
else
    echo "--- ss_clusters/ss_clusters.tsv already present, skipping SS clustering"
fi

# 4) Build the bilingual README_{EN,ES}.md from both result sets
echo "--- build README_EN.md + README_ES.md"
python3 "$REPO/scripts/build_cluster_readme.py" "$RUN_DIR" --name "$NAME"

echo "=== $(date -Is) done: $RUN_DIR"
