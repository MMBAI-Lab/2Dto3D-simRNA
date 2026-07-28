#!/usr/bin/env python3
"""
Render a PDB file as a simple 3D P-atom backbone trace PNG.

No PyMOL/ChimeraX available locally, so this uses Biopython + matplotlib's
mplot3d. The trace is drawn as a line connecting consecutive P atoms,
coloured 5'→3' with viridis. If P is missing (e.g. first residue), falls
back to C1'.

Usage: render_pdb_png.py <in.pdb> <out.png> [--title "..."] [--orient fixed|pca]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from matplotlib import cm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from Bio.PDB import PDBParser  # noqa: E402


BACKBONE_PRIORITY = ("P", "C1'", "C4'")


def backbone_coords(pdb_path: Path) -> np.ndarray:
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", str(pdb_path))
    coords = []
    for model in struct:
        for chain in model:
            for residue in chain:
                for name in BACKBONE_PRIORITY:
                    if residue.has_id(name):
                        coords.append(residue[name].coord)
                        break
        break  # first model only
    if not coords:
        raise RuntimeError(f"no backbone atoms found in {pdb_path}")
    return np.asarray(coords, dtype=float)


def orient_pca(xyz: np.ndarray) -> np.ndarray:
    xyz = xyz - xyz.mean(axis=0)
    _, _, vt = np.linalg.svd(xyz, full_matrices=False)
    return xyz @ vt.T


def render(pdb: Path, out: Path, title: str | None = None, orient: str = "pca") -> None:
    xyz = backbone_coords(pdb)
    if orient == "pca":
        xyz = orient_pca(xyz)
    else:
        xyz = xyz - xyz.mean(axis=0)

    n = len(xyz)
    colors = cm.viridis(np.linspace(0, 1, n))

    fig = plt.figure(figsize=(5, 5), dpi=150)
    ax = fig.add_subplot(111, projection="3d")

    # Segmented line coloured by residue index
    for i in range(n - 1):
        ax.plot(
            xyz[i:i + 2, 0], xyz[i:i + 2, 1], xyz[i:i + 2, 2],
            color=colors[i], linewidth=2.2,
        )
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=np.arange(n), cmap="viridis", s=14, depthshade=False)

    # Mark 5' and 3'
    ax.scatter(*xyz[0], color="#1b9e77", s=60, edgecolors="black", linewidths=0.8, label="5'")
    ax.scatter(*xyz[-1], color="#d95f02", s=60, edgecolors="black", linewidths=0.8, label="3'")
    ax.legend(loc="upper right", fontsize=8, frameon=False)

    # Equal aspect via manual bounds
    span = np.ptp(xyz, axis=0).max() / 2 + 2
    mid = xyz.mean(axis=0)
    ax.set_xlim(mid[0] - span, mid[0] + span)
    ax.set_ylim(mid[1] - span, mid[1] + span)
    ax.set_zlim(mid[2] - span, mid[2] + span)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.set_axis_off()
    ax.view_init(elev=20, azim=-60)

    if title:
        ax.set_title(title, fontsize=9)

    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pdb", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--title", default=None)
    ap.add_argument("--orient", choices=("fixed", "pca"), default="pca")
    args = ap.parse_args()
    render(args.pdb, args.out, args.title, args.orient)
    return 0


if __name__ == "__main__":
    sys.exit(main())
