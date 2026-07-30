#!/usr/bin/env python3
"""
Align and render every centroid PDB of a result set as a cartoon-ribbon PNG,
using ChimeraX headless.

Alignment is **per aptamer**, not per set. Least-squares C1' superposition is
only defined between runs of the same sequence, so centroids are grouped by the
sequence in `inputs/<set>/<run>/example`: each group gets its own anchor and its
own shared camera. Runs on different aptamers are never superposed onto each
other -- APT-PF1 happened to be a single-sequence set, which is why the earlier
version of this script could use one global reference.

Within a group the anchor is the most populated centroid (largest `pct` across
the group's `_centroids.tsv` files) -- the dominant conformer of that aptamer.

Effects per centroid:
  - `.pdb` is overwritten with the *aligned* coordinates (originals are
    regenerable from `rmsd_clusters/*.trafl`).
  - `.png` is rendered with the group's shared camera, so folds of the same
    aptamer can be compared without manually re-orienting.

Usage:
    python3 scripts/render_centroids_chimerax.py --set comercialApt
    python3 scripts/render_centroids_chimerax.py                # APT-PF1
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

CHIMERAX = "/home/pdans/chimerax/bin/chimerax-headless"

# Re-applied after each `open` so the newly loaded model gets the right style.
# (ChimeraX styling is per-model; nothing is inherited.)
STYLE = [
    "nucleotides #{m} tube/slab shape muffler",
    "cartoon style #{m} xsection round width 1.2 thickness 1.2",
    "rainbow #{m} residues palette rainbow",
]
SAVE_OPTS = "width 800 height 600 supersample 3 transparentBackground true"


def first_line(path: Path) -> str:
    if not path.exists():
        return ""
    for ln in path.read_text().splitlines():
        s = ln.strip()
        if s:
            return s
    return ""


def centroid_pct(pdb: Path) -> float:
    """Population of this centroid, read from its cutoff's _centroids.tsv."""
    tsv = pdb.parent / "_centroids.tsv"
    if not tsv.exists():
        return -1.0
    want = pdb.stem.replace("clust", "")
    with tsv.open() as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4 and parts[1] == want:
                return float(parts[3])
    return -1.0


def build_cxc(groups: list[tuple[str, Path, list[Path]]]) -> str:
    """One ChimeraX script covering every group; camera is re-derived per group."""
    lines: list[str] = ["set bgColor white", "lighting soft"]
    for label, anchor, others in groups:
        view_name = "cam_" + re.sub(r"[^A-Za-z0-9]", "_", label)
        lines.append(f"# ---- {label} ----")
        lines.append(f"open {anchor}")
        lines += [s.format(m=1) for s in STYLE]
        lines.append("view")
        lines.append(f"view name {view_name}")
        lines.append(f"save {anchor.with_suffix('.png')} {SAVE_OPTS}")

        for pdb in others:
            lines.append(f"open {pdb}")                     # -> model #2
            lines += [s.format(m=2) for s in STYLE]
            lines.append("align #2@C1' to #1@C1'")
            lines.append(f"save {pdb} models #2")           # aligned coords
            lines.append("hide #1 models")
            lines.append(f"view {view_name}")
            lines.append(f"save {pdb.with_suffix('.png')} {SAVE_OPTS}")
            lines.append("show #1 models")
            lines.append("close #2")
        lines.append("close #1")
    lines.append("exit")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--set", dest="set_name", default="APT-PF1")
    ap.add_argument("--chimerax", default=CHIMERAX)
    ap.add_argument("--ref-cutoff", action="append", default=[], metavar="RUN=CUTOFF",
                    help="Anchor each aptamer group on this run's dominant cluster at this "
                         "cutoff, instead of the most populated centroid overall. Repeatable. "
                         "Without it the anchor defaults to the largest pct, which always "
                         "lands on the coarsest cutoff and so is rarely a reported structure.")
    args = ap.parse_args()

    refcut: dict[str, float] = {}
    for spec in args.ref_cutoff:
        if "=" not in spec:
            print(f"--ref-cutoff expects RUN=CUTOFF, got {spec!r}", file=sys.stderr)
            return 1
        k, v = spec.split("=", 1)
        refcut[k] = float(v)

    global_root = args.repo / "results" / args.set_name / "global_analysis"
    in_root = args.repo / "inputs" / args.set_name
    res_root = args.repo / "results" / args.set_name
    if not global_root.is_dir():
        print(f"no centroids to render: {global_root} does not exist", file=sys.stderr)
        return 1
    if not Path(args.chimerax).exists():
        print(f"ChimeraX not found at {args.chimerax}", file=sys.stderr)
        return 1

    # Group centroid PDBs by the sequence of the run they came from.
    by_seq: dict[str, list[Path]] = {}
    seq_label: dict[str, str] = {}
    for run_dir in sorted(p for p in global_root.iterdir() if p.is_dir()):
        pdbs = sorted(run_dir.glob("cutoff_*A/clust*.pdb"))
        if not pdbs:
            continue
        seq = (first_line(in_root / run_dir.name / "example")
               or first_line(res_root / run_dir.name / "example")
               or run_dir.name)
        by_seq.setdefault(seq, []).extend(pdbs)
        seq_label.setdefault(seq, run_dir.name)

    if not by_seq:
        print(f"no centroid PDBs under {global_root}", file=sys.stderr)
        return 1

    def anchor_rank(pdb: Path) -> tuple[int, float]:
        """Prefer a centroid at its run's reference cutoff; then by population."""
        run = pdb.parent.parent.name
        want = refcut.get(run)
        m = re.match(r"cutoff_(\d+)A", pdb.parent.name)
        at_ref = 1 if (want is not None and m and float(m.group(1)) == want) else 0
        return (at_ref, centroid_pct(pdb))

    groups: list[tuple[str, Path, list[Path]]] = []
    total = 0
    for seq, pdbs in by_seq.items():
        anchor = max(pdbs, key=anchor_rank)
        others = [p for p in pdbs if p != anchor]
        label = seq_label[seq]
        groups.append((label, anchor, others))
        total += len(pdbs)
        print(f"group {label:<24} {len(seq):>3} nt  {len(pdbs):>3} centroids  "
              f"anchor={anchor.relative_to(global_root)} ({centroid_pct(anchor):.1f} %)")

    cxc_text = build_cxc(groups)
    with tempfile.NamedTemporaryFile("w", suffix=".cxc", delete=False) as f:
        f.write(cxc_text)
        cxc_path = f.name
    try:
        proc = subprocess.run([args.chimerax, cxc_path],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout[-2500:])
            sys.stderr.write(proc.stderr[-2500:])
            return proc.returncode
    finally:
        Path(cxc_path).unlink(missing_ok=True)

    missing = [p for grp in groups for p in ([grp[1]] + grp[2])
               if not p.with_suffix(".png").exists()]
    if missing:
        sys.stderr.write(f"missing PNGs: {len(missing)} of {total} (first: {missing[:3]})\n")
        return 1
    print(f"\nwrote {total} PNGs across {len(groups)} aptamer group(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
