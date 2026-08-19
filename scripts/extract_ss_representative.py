#!/usr/bin/env python3
"""
Extract a representative all-atom structure for a 2D (dot-bracket) pattern of a
run's low-T pool.

An RMSD cluster has a designated representative -- frame 1 of its `.trafl`. A 2D
pattern group has none: it is every frame whose detected base-pairing matches,
and those frames can differ widely in 3D. This picks the **lowest-energy** frame
of the group (field 4 of the `.trafl` header), which is reproducible and
defensible, but it is a representative of the *base-pairing*, not a centroid of
a 3D cluster. Do not read it as one.

Frame addressing: `ss_clusters/low_temp-NNNNNN.ss_detected` is the NNNNNN-th
frame of `low_temp.trafl` by position. The first field of a `.trafl` header is
the frame's original number in its source replica, not its position in the
pool, so only the positional index can be handed to SimRNA_trafl2pdbs.

Usage:
    python3 scripts/extract_ss_representative.py results/comercialApt/MC_J3T2_consensus \\
        --rank 1 --out /tmp/mc_rank1.pdb
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SS_FILE_RE = re.compile(r"low_temp-(\d+)\.ss_detected$")


def read_ss_table(run_dir: Path) -> list[tuple[int, str, int, float]]:
    tsv = run_dir / "ss_clusters" / "ss_clusters.tsv"
    rows = []
    with tsv.open() as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == 4:
                rows.append((int(p[0]), p[1], int(p[2]), float(p[3])))
    return rows


def frame_energies(trafl: Path) -> list[float]:
    """Energy per frame, indexed by position (0-based)."""
    out = []
    with trafl.open() as f:
        for line in f:
            if line and line[0].isdigit():
                parts = line.split()
                out.append(float(parts[3]) if len(parts) >= 4 else 0.0)
    return out


def frames_with_pattern(run_dir: Path, pattern: str) -> list[int]:
    """Positional (1-based) frame indices whose detected 2D equals `pattern`."""
    ss_dir = run_dir / "ss_clusters"
    hits = []
    with os.scandir(ss_dir) as it:
        for entry in it:
            m = SS_FILE_RE.search(entry.name)
            if not m:
                continue
            try:
                with open(entry.path) as fh:
                    for line in fh:
                        s = line.strip()
                        if s:
                            if s == pattern:
                                hits.append(int(m.group(1)))
                            break
            except OSError:
                continue
    return sorted(hits)


def extract_frame(run_dir: Path, index: int, out_pdb: Path, simrna_dir: Path) -> None:
    """Back-map one positional frame of low_temp.trafl to an all-atom PDB."""
    ref = run_dir / "ss_clusters" / "ref.pdb"
    trafl = run_dir / "low_temp.trafl"
    for p in (ref, trafl):
        if not p.exists():
            raise FileNotFoundError(p)

    env = os.environ.copy()
    env["PATH"] = f"{simrna_dir}:{env.get('PATH', '')}"
    with tempfile.TemporaryDirectory(prefix="ss_rep_") as td:
        tdp = Path(td)
        (tdp / "data").symlink_to(simrna_dir / "data")
        (tdp / "ref.pdb").symlink_to(ref.resolve())
        (tdp / "in.trafl").symlink_to(trafl.resolve())
        proc = subprocess.run(
            ["SimRNA_trafl2pdbs", "ref.pdb", "in.trafl", str(index), "AA"],
            cwd=td, env=env, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"SimRNA_trafl2pdbs failed for frame {index}\n"
                               f"{proc.stdout}\n{proc.stderr}")
        aa = sorted(tdp.glob("in-*_AA.pdb"))
        cg = [p for p in sorted(tdp.glob("in-*.pdb")) if not p.name.endswith("_AA.pdb")]
        produced = aa or cg
        if not produced:
            raise RuntimeError(f"no PDB produced for frame {index}: {list(tdp.iterdir())}")
        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(produced[0], out_pdb)


def pick(run_dir: Path, rank: int) -> tuple[str, float, int, float]:
    """Return (pattern, pct, best_frame_index, best_energy) for a 2D-pattern rank."""
    table = read_ss_table(run_dir)
    match = next((r for r in table if r[0] == rank), None)
    if match is None:
        raise SystemExit(f"rank {rank} not in {run_dir}/ss_clusters/ss_clusters.tsv")
    _, pattern, count, pct = match

    hits = frames_with_pattern(run_dir, pattern)
    if not hits:
        raise SystemExit(f"no frame carries pattern {pattern!r} in {run_dir}")
    if len(hits) != count:
        print(f"note: {len(hits)} frames match but the TSV counted {count}", file=sys.stderr)

    energies = frame_energies(run_dir / "low_temp.trafl")
    best = min(hits, key=lambda i: energies[i - 1] if i - 1 < len(energies) else 0.0)
    return pattern, pct, best, energies[best - 1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--rank", type=int, required=True, help="Rank in ss_clusters.tsv")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    args = ap.parse_args()

    pattern, pct, frame, energy = pick(args.run_dir, args.rank)
    print(f"rank {args.rank}: {pct:.2f} %  {pattern}")
    print(f"  lowest-energy frame: {frame}  (E = {energy:.3f})")
    extract_frame(args.run_dir, frame, args.out, args.repo / "SimRNA_64bitIntel_Linux")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
