#!/usr/bin/env python3
"""
Align every APT-PF1 centroid PDB to a single reference and re-render each as
a cartoon-ribbon PNG in that common coordinate frame.

The three APT-PF1 runs share the same 40-nt sequence, so per-C1' least-squares
alignment is well-defined. Choice of reference: the top-1 cluster @ 12 Å of
`DNA_as_RNA_VFold2D` (≈60 % of the low-T pool, reasonable anchor).

Effects per centroid:
  - `.pdb` is overwritten with the *aligned* coordinates (originals are
    regenerable from `rmsd_clusters/*.trafl`).
  - `.png` is re-rendered with a fixed camera shared across all 86 images, so
    folds can be visually compared without manually re-orienting.

Runs ChimeraX 1.11.1 headless via `/home/pdans/chimerax/bin/chimerax-headless`.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

CHIMERAX = "/home/pdans/chimerax/bin/chimerax-headless"
GLOBAL = Path("/home/pdans/Documents/2Dto3D-simRNA/results/APT-PF1/global_analysis")
REFERENCE = GLOBAL / "DNA_as_RNA_VFold2D" / "cutoff_12A" / "clust01.pdb"


# Commands re-applied after each `open` so the newly loaded model gets the
# right style. (ChimeraX styling is per-model; nothing is inherited.)
STYLE = [
    "nucleotides #{m} tube/slab shape muffler",
    "cartoon style #{m} xsection round width 1.2 thickness 1.2",
    "rainbow #{m} residues palette rainbow",
]


def build_cxc(reference: Path, others: list[Path]) -> str:
    lines: list[str] = []
    # --- Reference ---
    lines.append(f"open {reference}")
    lines += [s.format(m=1) for s in STYLE]
    lines.append("lighting soft")
    lines.append("set bgColor white")
    lines.append("view")
    lines.append("view name shared")  # save camera for reuse below

    ref_png = reference.with_suffix(".png")
    lines.append(f"save {ref_png} width 800 height 600 supersample 3 transparentBackground true")
    # (reference PDB coords are the anchor — no need to overwrite it.)

    # --- All other centroids ---
    for pdb in others:
        png = pdb.with_suffix(".png")
        lines.append(f"open {pdb}")                        # → model #2
        lines += [s.format(m=2) for s in STYLE]
        # Least-squares superposition on C1' atoms (shared between all 3 runs).
        lines.append("align #2@C1' to #1@C1'")
        # Save aligned coords (overwrites the original centroid PDB).
        lines.append(f"save {pdb} models #2")
        # Render with shared camera, reference hidden.
        lines.append("hide #1 models")
        lines.append("view shared")
        lines.append(f"save {png} width 800 height 600 supersample 3 transparentBackground true")
        lines.append("show #1 models")
        lines.append("close #2")

    lines.append("exit")
    return "\n".join(lines) + "\n"


def main() -> int:
    if not REFERENCE.exists():
        print(f"reference PDB missing: {REFERENCE}", file=sys.stderr)
        return 1
    pdbs = sorted(GLOBAL.glob("*/cutoff_*A/*.pdb"))
    others = [p for p in pdbs if p.resolve() != REFERENCE.resolve()]
    if not others:
        print("no other centroid PDBs found", file=sys.stderr)
        return 1
    print(f"aligning + rendering {len(others) + 1} centroids; ref = {REFERENCE.relative_to(GLOBAL)}")

    cxc_text = build_cxc(REFERENCE, others)
    with tempfile.NamedTemporaryFile("w", suffix=".cxc", delete=False) as f:
        f.write(cxc_text)
        cxc_path = f.name
    try:
        proc = subprocess.run(
            [CHIMERAX, cxc_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout[-2500:])
            sys.stderr.write(proc.stderr[-2500:])
            return proc.returncode
    finally:
        Path(cxc_path).unlink(missing_ok=True)

    missing = [p.with_suffix(".png") for p in pdbs if not p.with_suffix(".png").exists()]
    if missing:
        sys.stderr.write(f"missing PNGs: {len(missing)} (first: {missing[:3]})\n")
        return 1
    print(f"wrote {len(pdbs)} PNGs; aligned and overwrote {len(others)} PDBs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
