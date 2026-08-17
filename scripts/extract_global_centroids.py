#!/usr/bin/env python3
"""
For each run × cutoff × top-N cluster of a result set, extract the centroid
frame as an all-atom PDB (back-mapped by `SimRNA_trafl2pdbs`), also save its
detected `.ss_detected`, and render a PNG snapshot.

Defaults to the APT-PF1 set; pass `--set <name>` for any other, e.g.

    python3 scripts/extract_global_centroids.py --set comercialApt

Layout written under results/<set>/global_analysis/<run>/cutoff_<N>A/:
  clustNN.pdb          — AA centroid (frame 1 of cluster .trafl)
  clustNN.ss_detected  — detected dot-bracket
  clustNN.png          — rendered backbone trace
  _centroids.tsv       — (rank, cutoff, cluster_id, frames, pct, ss)

Expected inputs (per run_dir):
  <run_dir>/<name>_01-000001.pdb  (reference for trafl2pdbs)
  <run_dir>/rmsd_clusters/low_temp_thrs<X>A_clust<NN>.trafl
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

CLUSTER_RE = re.compile(r"low_temp_thrs(?P<cutoff>[\d.]+)A_clust(?P<clust>\d+)\.trafl")


def count_frames(trafl: Path) -> int:
    n = 0
    with trafl.open() as f:
        for line in f:
            if line and line[0].isdigit():
                n += 1
    return n


def read_ss(ss_path: Path) -> str:
    if not ss_path.exists():
        return ""
    for line in ss_path.read_text().splitlines():
        s = line.strip()
        if s and all(c in ".()[]{}<>-&" for c in s):
            return s
    return ""


def discover(rmsd_dir: Path) -> dict[float, list[tuple[int, int, Path]]]:
    """Return {cutoff -> [(cluster_id, frame_count, trafl_path)] desc by count}."""
    by_cut: dict[float, list[tuple[int, int, Path]]] = {}
    for path in rmsd_dir.glob("low_temp_thrs*A_clust*.trafl"):
        m = CLUSTER_RE.match(path.name)
        if not m:
            continue
        cutoff = float(m.group("cutoff"))
        clust_id = int(m.group("clust"))
        by_cut.setdefault(cutoff, []).append((clust_id, count_frames(path), path))
    return {c: sorted(lst, key=lambda t: (-t[1], t[0])) for c, lst in by_cut.items()}


def backmap_frame1(
    ref_pdb: Path,
    cluster_trafl: Path,
    out_pdb: Path,
    out_ss: Path,
    simrna_dir: Path,
    render_script: Path,
    png_path: Path,
    title: str,
    render: bool = True,
) -> None:
    """Run SimRNA_trafl2pdbs (AA) for frame 1 in a tempdir, move outputs, render PNG.

    Back-mapping needs only the SimRNA binary. Rendering needs numpy/matplotlib/
    biopython, which is why it is separable: `render=False` still produces every
    centroid PDB, and scripts/render_centroids_chimerax.py can draw them later
    with no Python dependencies at all.
    """
    env = os.environ.copy()
    env["PATH"] = f"{simrna_dir}:{env.get('PATH', '')}"
    with tempfile.TemporaryDirectory(prefix="trafl2pdbs_") as td:
        tdp = Path(td)
        (tdp / "data").symlink_to(simrna_dir / "data")
        ref_link = tdp / "ref.pdb"
        trafl_link = tdp / "in.trafl"
        ref_link.symlink_to(ref_pdb.resolve())
        trafl_link.symlink_to(cluster_trafl.resolve())

        # SimRNA_trafl2pdbs ref in.trafl 1 AA
        proc = subprocess.run(
            ["SimRNA_trafl2pdbs", "ref.pdb", "in.trafl", "1", "AA"],
            cwd=td, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"SimRNA_trafl2pdbs failed for {cluster_trafl}\n"
                f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
            )

        aa_candidates = sorted(tdp.glob("in-*_AA.pdb"))
        cg_candidates = [p for p in sorted(tdp.glob("in-*.pdb")) if not p.name.endswith("_AA.pdb")]
        produced_pdb = aa_candidates or cg_candidates
        if not produced_pdb:
            raise RuntimeError(f"no PDB produced for {cluster_trafl}; tmp contents: {list(tdp.iterdir())}")
        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(produced_pdb[0], out_pdb)
        produced_ss = sorted(tdp.glob("in-*.ss_detected"))
        if produced_ss:
            shutil.copy(produced_ss[0], out_ss)

    if not render:
        return

    # Render PNG
    render_proc = subprocess.run(
        [sys.executable, str(render_script), str(out_pdb), str(png_path), "--title", title],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    if render_proc.returncode != 0:
        sys.stderr.write(f"render failed for {out_pdb}:\n{render_proc.stderr}\n")


def process_run(
    run_dir: Path,
    global_dir: Path,
    simrna_dir: Path,
    render_script: Path,
    name: str,
    top_n: int,
    render: bool = True,
) -> None:
    rmsd_dir = run_dir / "rmsd_clusters"
    ref_pdb = run_dir / f"{name}_01-000001.pdb"
    if not ref_pdb.exists():
        raise FileNotFoundError(f"missing reference PDB: {ref_pdb}")
    out_run = global_dir / run_dir.name
    out_run.mkdir(parents=True, exist_ok=True)

    by_cutoff = discover(rmsd_dir)

    for cutoff, clusters in sorted(by_cutoff.items()):
        total = sum(c for _, c, _ in clusters) or 1
        top = clusters[:top_n]
        cut_dir = out_run / f"cutoff_{int(round(cutoff))}A"
        cut_dir.mkdir(parents=True, exist_ok=True)
        tsv = cut_dir / "_centroids.tsv"
        with tsv.open("w") as tf:
            tf.write("rank\tcluster\tframes\tpct\tss\n")
            for rank, (cid, count, trafl) in enumerate(top, start=1):
                pct = 100.0 * count / total
                pdb_out = cut_dir / f"clust{cid:02d}.pdb"
                ss_out = cut_dir / f"clust{cid:02d}.ss_detected"
                png_out = cut_dir / f"clust{cid:02d}.png"
                title = f"{run_dir.name}  cutoff {cutoff:.0f}Å  rank {rank} (clust {cid:02d}, {pct:.1f}%)"
                print(f"[{run_dir.name}] cutoff={cutoff:.0f}Å rank={rank} clust={cid:02d} ({count} frames, {pct:.2f}%)")
                if not pdb_out.exists():
                    backmap_frame1(
                        ref_pdb, trafl, pdb_out, ss_out, simrna_dir, render_script, png_out,
                        title, render,
                    )
                elif render and not png_out.exists():
                    subprocess.run(
                        [sys.executable, str(render_script), str(pdb_out), str(png_out), "--title", title],
                        check=False,
                    )
                ss = read_ss(ss_out)
                tf.write(f"{rank}\t{cid:02d}\t{count}\t{pct:.2f}\t{ss}\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--set", dest="set_name", default="APT-PF1",
                    help="Result set under results/ (default: APT-PF1)")
    ap.add_argument("--runs", nargs="*", default=None,
                    help="Run directory names to process. Default: every subdirectory of "
                         "results/<set>/ that has a rmsd_clusters/ directory.")
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--name", default="example")
    ap.add_argument("--no-render", action="store_true",
                    help="Only back-map the centroid PDBs; skip the matplotlib PNG pass. "
                         "Use when numpy/matplotlib/biopython are unavailable, or when "
                         "render_centroids_chimerax.py will draw them instead.")
    args = ap.parse_args()

    simrna_dir = args.repo / "SimRNA_64bitIntel_Linux"
    render_script = args.repo / "scripts" / "render_pdb_png.py"
    set_root = args.repo / "results" / args.set_name
    global_dir = set_root / "global_analysis"

    if not set_root.is_dir():
        print(f"no such result set: {set_root}", file=sys.stderr)
        return 1

    if args.runs:
        runs = list(args.runs)
    else:
        runs = sorted(
            p.name for p in set_root.iterdir()
            if p.is_dir() and (p / "rmsd_clusters").is_dir()
        )
    if not runs:
        print(f"no runs with rmsd_clusters/ found under {set_root}", file=sys.stderr)
        return 1

    global_dir.mkdir(parents=True, exist_ok=True)
    print(f"set '{args.set_name}': {len(runs)} run(s) -> {', '.join(runs)}")
    for r in runs:
        process_run(set_root / r, global_dir, simrna_dir, render_script, args.name,
                    args.top_n, render=not args.no_render)

    print(f"\nCentroids written under {global_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
