#!/usr/bin/env python3
"""
Compute per-frame distances between two C1' atoms from a SimRNA trajectory.

This script mirrors the CLI/UX conventions of `sasa_from_simrna.py`, adapted for
C1'–C1' distances. It supports two modes:

1) All-atom (AA) mode (preferred):
   - Uses `SimRNA_trafl2pdbs` to reconstruct AA PDBs for selected frames.
   - Parses each frame PDB, finds C1' (or rare alias C1*) for two selected residues,
     computes Euclidean distance per frame, and writes results.

2) Coarse-grained (CG) fallback mode:
   - Parses `.trafl` directly. If `--c4prime-proxy` is set, measures distance between
     the C4' beads (bead #2 per 5-bead residue block; index 1) of the two residues
     as an approximation of the C1' distance. Otherwise, errors with guidance to
     use `--aa` or `--c4prime-proxy`.

Outputs:
  * Per-frame CSV (default: c1prime_distance_per_frame.csv)
  * Optional run-row appends (wide-run table)
  * Optional plot (distance vs. frame + boxplot summary)
  * Summary CSV with descriptive statistics (median/mean/std/min/max/q1/q3)


Examples
========
AA mode for all frames:

    python distance.py \
        --pdb reference.pdb --trafl traj.trafl \
        --nt1 A:10 --nt2 B:25 \
        --aa --frames ":" \
        --out per_frame.csv --plot distances.png

CG proxy (if AA unavailable):

    python distance.py \
        --pdb reference.pdb --trafl traj.trafl \
        --nt1 10 --nt2 25 \
        --c4prime-proxy \
        --out per_frame.csv --plot distances.png

Testing notes (see bottom of file for mini tests and guidelines):
  - Unit-test the PDB parser on tiny snippets containing C1' and C1*.
  - Test selector parsing ("A:10", "10") and ambiguous/no-chain cases.
  - Test CG indexing assumes bead order [P, C4', N, B1, B2]; use bead index 1.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from tqdm import tqdm  # type: ignore
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    def tqdm(x: Iterable, **_: dict):
        return x

# ------------------------------
# Helpers for CLI ergonomics
# ------------------------------

Selector = Tuple[Optional[str], int]  # (chain, resid)


@dataclass
class Args:
    pdb: Path
    trafl: Path
    nt1: str
    nt2: str
    aa: bool
    c4prime_proxy: bool
    simrna_bin: Optional[Path]
    frames: str
    out: Path
    plot: Optional[Path]
    append_run_rows: Optional[Path]
    boxplot_csv: Optional[Path]
    run_id: Optional[str]
    keep_pdb: bool


CHECK = "\N{WHITE HEAVY CHECK MARK}"
CROSS = "\N{CROSS MARK}"
INFO = "\N{INFORMATION SOURCE}"


# ------------------------------
# SimRNA trafl2pdbs detection and runner (from sasa.py)
# ------------------------------

def detect_simrna_trafl2pdbs(provided: Optional[str] = None) -> Optional[str]:
    """Detect SimRNA_trafl2pdbs executable. Return path or None.
    By default we look on PATH (user requested default behaviour).
    """
    if provided:
        if os.path.isfile(provided) and os.access(provided, os.X_OK):
            return provided
        return None
    name = shutil.which('SimRNA_trafl2pdbs')
    return name


def run_simrna_trafl2pdbs(simrna_bin: str, reference_pdb: str, trafl: str, frames: str, aa: bool, out_dir: str) -> List[str]:
    """Call SimRNA_trafl2pdbs to convert desired frames into PDB files.
    Returns list of generated pdb file paths (in increasing frame order found).
    """
    os.makedirs(out_dir, exist_ok=True)
    cmd = [simrna_bin, reference_pdb, trafl, frames]
    if aa:
        cmd.append('AA')
    print(f"Running: {' '.join(cmd)}")
    
    # Run in current directory (where data symlink exists)
    subprocess.check_call(cmd)
    
    # Look for generated files in current directory
    base = os.path.splitext(os.path.basename(trafl))[0]
    reference_pdb_name = os.path.basename(reference_pdb)
    
    # Find all generated files (but not the reference PDB)
    all_generated_files = [f for f in os.listdir('.') if f.startswith(base) and f != reference_pdb_name]
    
    # Separate AA PDBs, CG PDBs, and ss_detected files
    aa_pdbs = []
    cg_pdbs = []
    ss_detected_files = []
    
    for f in all_generated_files:
        if f.endswith('_AA.pdb'):
            aa_pdbs.append(f)
        elif f.endswith('.pdb') and not f.endswith('_AA.pdb'):
            # This is a CG PDB (frame output without _AA suffix)
            if '-' in f:  # Make sure it's a frame file, not the reference
                cg_pdbs.append(f)
        elif f.endswith('.ss_detected'):
            ss_detected_files.append(f)
    
    files_to_process = []
    
    if aa:
        # For AA mode, move only the AA PDBs to temp directory
        for f in sorted(aa_pdbs):
            src = f
            dst = os.path.join(out_dir, f)
            try:
                shutil.move(src, dst)
                files_to_process.append(dst)
            except Exception as e:
                print(f"Warning: Could not move {src} to {dst}: {e}")
        
        # Clean up CG PDBs (we don't need them in AA mode)
        for f in cg_pdbs:
            try:
                os.remove(f)
                print(f"Cleaned up CG PDB: {f}")
            except Exception as e:
                print(f"Warning: Could not remove CG PDB {f}: {e}")
    else:
        # For CG mode, move CG PDBs to temp directory
        for f in sorted(cg_pdbs):
            src = f
            dst = os.path.join(out_dir, f)
            try:
                shutil.move(src, dst)
                files_to_process.append(dst)
            except Exception as e:
                print(f"Warning: Could not move {src} to {dst}: {e}")
    
    # Always clean up .ss_detected files (we don't use them in SASA calculation)
    for f in ss_detected_files:
        try:
            os.remove(f)
            print(f"Cleaned up .ss_detected file: {f}")
        except Exception as e:
            print(f"Warning: Could not remove .ss_detected file {f}: {e}")
    
    print(f"Found and moved {len(files_to_process)} {'AA' if aa else 'CG'} PDB files to temporary directory")
    return files_to_process


# ------------------------------
# Parsing helpers
# ------------------------------

def parse_selector(sel: str) -> Selector:
    """Parse residue selector like 'A:10' or '10'. Returns (chain or None, resid int).
    Chain may be empty string; treat as None.
    """
    m = re.match(r"^(?:(?P<chain>[^:]+):)?(?P<resid>\d+)$", sel.strip())
    if not m:
        raise ValueError(f"Invalid selector '{sel}'. Use A:10 or 10.")
    chain = m.group("chain")
    if chain == "" or chain is None:
        chain = None
    resid = int(m.group("resid"))
    return (chain, resid)


def parse_frames_arg(frames: str, total_frames_hint: Optional[int] = None) -> List[int]:
    """Parse --frames syntax similar to reference.

    Accepted forms:
      ':'              -> all frames (use total_frames_hint if provided, else empty -> later treat as all)
      '5'              -> single frame 5
      '1:10'           -> inclusive range [1,10]
      '1:10:2'         -> step (1,3,5,...)
      '1,5,8'          -> explicit list (commas)
    Returns a sorted list of 1-based frame indices (unique).
    """
    s = frames.strip()
    if s == ":":
        if total_frames_hint is None:
            return []  # empty means "all"; caller may resolve later
        return list(range(1, total_frames_hint + 1))
    if "," in s:
        vals = sorted({int(x.strip()) for x in s.split(",") if x.strip()})
        return vals
    if ":" in s:
        parts = [p for p in s.split(":")]
        if len(parts) == 2:
            a, b = int(parts[0]), int(parts[1])
            return list(range(a, b + 1))
        elif len(parts) == 3:
            a, b, c = int(parts[0]), int(parts[1]), int(parts[2])
            return list(range(a, b + 1, c))
        else:
            raise ValueError(f"Invalid --frames format: {frames}")
    # single integer
    return [int(s)]


AtomKey = Tuple[str, int]  # (chain, resid)


def parse_pdb_find_c1prime(pdb_path: Path) -> Dict[AtomKey, np.ndarray]:
    """Return mapping (chain,resid) -> C1' 3D coordinate for a PDB frame.

    Recognizes atom names "C1'" and rare alias "C1*" (maps to same). Only ATOM/HETATM lines.
    """
    c1_map: Dict[AtomKey, np.ndarray] = {}
    with pdb_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            # PDB fixed columns (1-based):
            #  1-6 rec, 7-11 serial, 13-16 atom name, 17 altLoc, 18-20 resName,
            # 22 chainID, 23-26 resSeq, 31-38 x, 39-46 y, 47-54 z
            atom_name = line[12:16].strip()
            if atom_name in {"C1'", "C1*"}:
                chain = line[21].strip() or ""
                try:
                    resid = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                key: AtomKey = (chain, resid)
                c1_map[key] = np.array([x, y, z], dtype=float)
    return c1_map


# ------------------------------
# TRAFL parsing (CG coords)
# ------------------------------

def parse_trafl_frames(trafl: Path) -> Tuple[np.ndarray, int, int]:
    """Parse a .trafl file into an array of shape (n_frames, n_atoms, 3).

    Follows the common SimRNA format: pairs of lines per frame (header + coords).
    Returns (coords, n_frames, n_atoms).
    """
    lines = trafl.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) % 2 != 0:
        raise RuntimeError(".trafl file malformed: expected pairs of lines per frame.")
    frames = len(lines) // 2
    coords_list: List[np.ndarray] = []
    n_atoms: Optional[int] = None
    for i in range(frames):
        head = lines[2 * i].strip()
        xyzs = lines[2 * i + 1].strip().split()
        if len(xyzs) % 3 != 0:
            raise RuntimeError(f"Frame {i+1}: coordinate line length not divisible by 3.")
        fa = np.array(list(map(float, xyzs)), dtype=float).reshape(-1, 3)
        if n_atoms is None:
            n_atoms = fa.shape[0]
        elif fa.shape[0] != n_atoms:
            raise RuntimeError("Inconsistent atom counts between frames in .trafl")
        coords_list.append(fa)
    assert n_atoms is not None
    coords = np.stack(coords_list, axis=0)
    return coords, frames, n_atoms


# ------------------------------
# Core analysis
# ------------------------------

def _resolve_atom_for_selector(
    c1_map: Dict[AtomKey, np.ndarray], sel: Selector
) -> Tuple[str, int, np.ndarray]:
    """Resolve a selector (chain,resid) against c1_map.

    If chain is None, match first residue with resid across chains; warn if multiple chains.
    Returns canonical (chain,resid, coord). Chain is a single-character string (may be '').
    """
    chain, resid = sel
    matches: List[Tuple[str, int, np.ndarray]] = []
    if chain is None:
        for (ch, rs), coord in c1_map.items():
            if rs == resid:
                matches.append((ch, rs, coord))
        if not matches:
            raise KeyError(f"Resid {resid} not found in any chain for this frame.")
        if len(matches) > 1:
            chosen = sorted(matches, key=lambda t: (t[0], t[1]))[0]
            print(
                f"{INFO} Ambiguous resid {resid}: multiple chains present; using first match '{chosen[0] or ':'}:{chosen[1]}'"
            )
            return chosen
        return matches[0]
    else:
        key = (chain, resid)
        if key not in c1_map:
            # Show available chain/resid combos for debugging
            raise KeyError(
                f"Resid {chain}:{resid} not found (no C1' atom) in this frame."
            )
        return (chain, resid, c1_map[key])


def run_analysis(args: Args) -> pd.DataFrame:
    sel1 = parse_selector(args.nt1)
    sel2 = parse_selector(args.nt2)

    mode_used: Optional[str] = None

    if args.aa:
        simrna_bin = detect_simrna_trafl2pdbs(str(args.simrna_bin) if args.simrna_bin else None)
        if not simrna_bin:
            print(
                f"{INFO} SimRNA_trafl2pdbs not found; falling back to CG mode. Try setting --simrna-bin or add SimRNA to PATH."
            )
        else:
            # Use the robust run_simrna_trafl2pdbs function from sasa.py
            print("Converting trajectory frames to all-atom PDBs...")
            tempdir = tempfile.mkdtemp(prefix='simrna_pdbs_distance_')
            frames_spec = args.frames if args.frames else ':'
            
            try:
                pdb_frames = run_simrna_trafl2pdbs(
                    simrna_bin, str(args.pdb), str(args.trafl), frames_spec, aa=True, out_dir=tempdir
                )
                if len(pdb_frames) == 0:
                    print('No PDB frames produced by converter; falling back to CG mode')
                    args.aa = False
                else:
                    print(f"✓ Converted {len(pdb_frames)} frames to all-atom PDBs")
                    mode_used = "aa"
                    
                    # Process frames with progress bar
                    rows: List[Dict[str, object]] = []
                    frame_iterator = tqdm(enumerate(sorted(pdb_frames)), total=len(pdb_frames), desc="Processing frames") if HAS_TQDM else enumerate(sorted(pdb_frames))
                    
                    for idx, pdb_path in frame_iterator:
                        c1_map = parse_pdb_find_c1prime(Path(pdb_path))
                        if not c1_map:
                            raise RuntimeError(
                                f"No C1' atoms found in {os.path.basename(pdb_path)}; ensure AA reconstruction produced sugar atoms."
                            )
                        ch1, rs1, v1 = _resolve_atom_for_selector(c1_map, sel1)
                        ch2, rs2, v2 = _resolve_atom_for_selector(c1_map, sel2)
                        dist = float(np.linalg.norm(v1 - v2))
                        rows.append(
                            {
                                "frame": idx + 1,
                                "nt1": f"{ch1 or ''}:{rs1}" if ch1 != "" else f":{rs1}",
                                "nt2": f"{ch2 or ''}:{rs2}" if ch2 != "" else f":{rs2}",
                                "distance_A": dist,
                                "mode": mode_used,
                            }
                        )
                        
                        # Remove PDB if not keeping
                        if not args.keep_pdb:
                            try:
                                os.remove(pdb_path)
                            except Exception:
                                pass
                    
                    df = pd.DataFrame(rows)
                    print(f"{CHECK} Computed distances for {len(df)} frames in AA mode.")
                    
                    # Clean up tempdir if not keeping
                    if not args.keep_pdb:
                        try:
                            shutil.rmtree(tempdir, ignore_errors=True)
                        except Exception:
                            pass
                    
                    return df
                    
            except Exception as e:
                print(f"{INFO} AA path failed ({e}). Falling back to CG path if possible...")
                # Clean up on failure
                try:
                    if not args.keep_pdb:
                        shutil.rmtree(tempdir, ignore_errors=True)
                except Exception:
                    pass

    # CG path
    coords, n_frames, n_atoms = parse_trafl_frames(args.trafl)
    if args.frames and args.frames.strip() != ":":
        frames_sel = parse_frames_arg(args.frames, total_frames_hint=n_frames)
        # Convert to 0-based indices
        frame_indices = [i - 1 for i in frames_sel]
        coords = coords[frame_indices]
        n_frames = coords.shape[0]
    elif args.frames and args.frames.strip() == ":":
        # all frames already loaded
        pass

    # Validate 5 beads per residue
    if n_atoms % 5 != 0:
        raise RuntimeError(
            f".trafl appears not to be 5-bead-per-residue (atoms/frame={n_atoms})."
        )

    if not args.c4prime_proxy:
        raise RuntimeError(
            "C1' atoms exist only in AA models. Rerun with --aa to reconstruct all-atom PDBs, "
            "or use --c4prime-proxy to approximate using C4' beads in CG mode."
        )

    n_res = n_atoms // 5
    # Bead order per residue: [P, C4', N, B1, B2]; choose bead index 1 for C4'
    # Residues are 1-based in selectors; map resid r -> bead index base = (r-1)*5

    def bead_index_for_residue(resid_1b: int) -> int:
        return (resid_1b - 1) * 5 + 1  # 0-based atom index for C4'

    (ch1, resid1) = sel1
    (ch2, resid2) = sel2
    if ch1 is not None or ch2 is not None:
        print(
            f"{INFO} Chains are ignored in CG mode; using resid indices only (resid1={resid1}, resid2={resid2})."
        )

    if not (1 <= resid1 <= n_res) or not (1 <= resid2 <= n_res):
        raise RuntimeError(
            f"Resid indices out of range for CG beads (n_res={n_res}, requested {resid1} and {resid2})."
        )

    idx1 = bead_index_for_residue(resid1)
    idx2 = bead_index_for_residue(resid2)

    dists = np.linalg.norm(coords[:, idx1, :] - coords[:, idx2, :], axis=1)

    # Compose DataFrame
    df = pd.DataFrame(
        {
            "frame": np.arange(1, n_frames + 1, dtype=int),
            "nt1": [f":{resid1}"] * n_frames,
            "nt2": [f":{resid2}"] * n_frames,
            "distance_A": dists.astype(float),
            "mode": ["cg_c4prime_proxy"] * n_frames,
        }
    )
    print(
        f"{CHECK} Computed distances for {n_frames} frames in CG C4' proxy mode (approximation)."
    )
    return df


# ------------------------------
# I/O: CSVs and plots
# ------------------------------

def write_per_frame_csv(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"{CHECK} Wrote per-frame CSV: {out_path}")


def append_run_rows(df: pd.DataFrame, path: Path, run_id: Optional[str]) -> None:
    if run_id is None:
        run_id = "run"
    rows = df[["frame", "nt1", "nt2", "distance_A", "mode"]].copy()
    rows.insert(0, "run_id", run_id)
    header = not path.exists()
    with path.open("a", newline="") as fh:
        writer = csv.writer(fh)
        if header:
            writer.writerow(["run_id", "mode", "frame", "nt1", "nt2", "distance_A"])
        for _, r in rows.iterrows():
            writer.writerow([r["run_id"], r["mode"], int(r["frame"]), r["nt1"], r["nt2"], float(r["distance_A"])])
    print(f"{CHECK} Appended run rows to: {path}")


def summarize_and_write_boxplot_csv(
    df: pd.DataFrame, path: Path, run_id: Optional[str]
) -> pd.DataFrame:
    vals = df["distance_A"].astype(float).to_numpy()
    if vals.size == 0:
        raise RuntimeError("No values to summarize.")
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    summary = {
        "run_id": run_id or "run",
        "nt1": df.iloc[0]["nt1"],
        "nt2": df.iloc[0]["nt2"],
        "median": float(med),
        "mean": float(np.mean(vals)),
        "std": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
        "min": float(np.min(vals)),
        "max": float(np.max(vals)),
        "q1": float(q1),
        "q3": float(q3),
    }
    header = not path.exists()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "run_id",
                "nt1",
                "nt2",
                "median",
                "mean",
                "std",
                "min",
                "max",
                "q1",
                "q3",
            ],
        )
        if header:
            writer.writeheader()
        writer.writerow(summary)
    print(f"{CHECK} Wrote/updated boxplot summary CSV: {path}")
    return pd.DataFrame([summary])


def make_plot(df: pd.DataFrame, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(12, 6))

    # Panel 1: line plot (distance vs frame)
    ax1 = fig.add_subplot(1, 2, 1)
    
    # Determine alpha and marker size based on number of points
    n_points = len(df)
    if n_points > 10000:
        alpha = 0.3
        markersize = 1
        linewidth = 0.5
    elif n_points > 5000:
        alpha = 0.4
        markersize = 1.5
        linewidth = 0.7
    elif n_points > 1000:
        alpha = 0.6
        markersize = 2
        linewidth = 0.8
    else:
        alpha = 0.8
        markersize = 3
        linewidth = 1
    
    # Plot with enhanced visibility
    ax1.plot(df["frame"].to_numpy(), df["distance_A"].to_numpy(), 
             marker="o", markersize=markersize, alpha=alpha, linewidth=linewidth,
             color='steelblue', markeredgewidth=0)
    
    # Add a smoothed trend line for dense data
    if n_points > 500:
        from scipy import ndimage
        smoothed = ndimage.gaussian_filter1d(df["distance_A"].to_numpy(), sigma=max(1, n_points//200))
        ax1.plot(df["frame"].to_numpy(), smoothed, 
                 color='red', linewidth=2, alpha=0.8, label='Smoothed trend')
        ax1.legend()
    
    ax1.set_xlabel("Frame (1-based)")
    ax1.set_ylabel("C1'–C1' Distance (Å)")
    ax1.set_title(f"Per-frame distance (n={n_points})")
    ax1.grid(True, alpha=0.3)
    
    # Add statistics text box
    dist_values = df["distance_A"].to_numpy()
    stats_text = f'Mean: {np.mean(dist_values):.2f} Å\nStd: {np.std(dist_values):.2f} Å\nRange: {np.min(dist_values):.2f}-{np.max(dist_values):.2f} Å'
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Panel 2: enhanced boxplot with violin plot overlay
    ax2 = fig.add_subplot(1, 2, 2)
    
    # Create violin plot for better distribution visualization
    parts = ax2.violinplot([dist_values], positions=[1], widths=[0.6], 
                          showmeans=True, showmedians=True, showextrema=False)
    
    # Style the violin plot
    for pc in parts['bodies']:
        pc.set_facecolor('lightblue')
        pc.set_alpha(0.6)
    
    # Overlay boxplot
    bp = ax2.boxplot([dist_values], positions=[1], widths=[0.4], 
                     patch_artist=True, notch=True)
    
    # Style the boxplot
    bp['boxes'][0].set_facecolor('steelblue')
    bp['boxes'][0].set_alpha(0.8)
    
    # Add scatter points for outliers if not too many
    if n_points <= 1000:
        # Add jittered points
        y_vals = dist_values
        x_vals = np.random.normal(1, 0.04, size=len(y_vals))
        ax2.scatter(x_vals, y_vals, alpha=min(0.6, 500/n_points), s=8, color='darkblue')
    
    ax2.set_xlim(0.5, 1.5)
    ax2.set_xticks([1])
    ax2.set_xticklabels(['Distance'])
    ax2.set_ylabel("C1'–C1' Distance (Å)")
    ax2.set_title("Distribution summary")
    ax2.grid(True, alpha=0.3, axis='y')

    # Overall title with mode information
    mode = df.iloc[0]['mode'] if 'mode' in df.columns else 'unknown'
    nt1 = df.iloc[0]['nt1'] if 'nt1' in df.columns else 'NT1'
    nt2 = df.iloc[0]['nt2'] if 'nt2' in df.columns else 'NT2'
    
    fig.suptitle(f"C1' distance: {nt1} ↔ {nt2} ({mode} mode)", fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    # Save with high DPI for better quality
    fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"{CHECK} Saved enhanced plot: {out_path}")


# ------------------------------
# CLI
# ------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compute per-frame distances between two C1' atoms from a SimRNA trajectory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--pdb", required=True, type=Path, help="Reference first-frame PDB (AA PDB)")
    p.add_argument("--trafl", required=True, type=Path, help="Path to .trafl trajectory")
    p.add_argument("--nt1", required=True, help="Residue selector for nucleotide 1 (A:10 or 10)")
    p.add_argument("--nt2", required=True, help="Residue selector for nucleotide 2 (A:10 or 10)")

    p.add_argument("--aa", action="store_true", help="Reconstruct AA PDBs per frame using SimRNA_trafl2pdbs")
    p.add_argument(
        "--c4prime-proxy",
        action="store_true",
        help="In CG mode, measure distance using C4' bead (bead #2 per residue) as proxy",
    )
    p.add_argument(
        "--simrna-bin",
        type=Path,
        default=None,
        help="Path to SimRNA_trafl2pdbs (overrides PATH lookup)",
    )
    p.add_argument(
        "--frames",
        default=":",
        help="Frames to process (':', 'N', 'A:B', 'A:B:S', or 'i,j,k')",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("c1prime_distance_per_frame.csv"),
        help="Per-frame CSV output path",
    )
    p.add_argument("--plot", type=Path, default=None, help="Optional plot PNG path")
    p.add_argument(
        "--append-run-rows",
        type=Path,
        default=None,
        help="Append wide run rows to this CSV (columns: run_id, mode, frame, nt1, nt2, distance_A)",
    )
    p.add_argument(
        "--boxplot-csv",
        type=Path,
        default=None,
        help="Summary CSV path (median/mean/std/min/max/q1/q3); defaults to <out>.boxplot_summary.csv",
    )
    p.add_argument("--run-id", default=None, help="Run identifier to include in appended rows/summary")
    p.add_argument(
        "--keep-pdb",
        action="store_true",
        help="Keep generated AA PDBs (default: delete tempdir)",
    )
    return p


# ------------------------------
# Main
# ------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    ns = parser.parse_args(argv)
    args = Args(
        pdb=ns.pdb,
        trafl=ns.trafl,
        nt1=ns.nt1,
        nt2=ns.nt2,
        aa=ns.aa,
        c4prime_proxy=ns.c4prime_proxy,
        simrna_bin=ns.simrna_bin,
        frames=ns.frames,
        out=ns.out,
        plot=ns.plot,
        append_run_rows=ns.append_run_rows,
        boxplot_csv=ns.boxplot_csv,
        run_id=ns.run_id,
        keep_pdb=ns.keep_pdb,
    )

    # Basic validations
    if not args.pdb.exists():
        print(f"{CROSS} Missing --pdb file: {args.pdb}")
        return 1
    if not args.trafl.exists():
        print(f"{CROSS} Missing --trafl file: {args.trafl}")
        return 1

    try:
        df = run_analysis(args)
        write_per_frame_csv(df, args.out)
        if args.append_run_rows:
            append_run_rows(df, args.append_run_rows, run_id=args.run_id)
        # Summary CSV
        box_path = args.boxplot_csv or Path(str(args.out) + ".boxplot_summary.csv")
        summarize_and_write_boxplot_csv(df, box_path, run_id=args.run_id)
        # Plot
        if args.plot:
            make_plot(df, args.plot)
        print(f"{CHECK} Done.")
        return 0
    except Exception as e:
        print(f"{CROSS} Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())

# ------------------------------
# Mini tests & guidance (comments only)
# ------------------------------
# - Test PDB parser with lines like:
#   ATOM      1  C1'  A  A  10      10.000  11.000  12.000  1.00 20.00           C
#   ATOM      2  C1*  G  B  25      20.000  21.000  22.000  1.00 20.00           C
#   Ensure both map to keys ('A',10) and ('B',25) with correct coords.
# - Test parse_selector: parse_selector('A:10')->('A',10); parse_selector('10')->(None,10)
# - Test CG indexing: for n_res=3, resid=2 -> base=(2-1)*5=5; C4' index=6th atom 0-based -> 5+1
# - Test --frames parsing variants: ':', '5', '1:3', '1:5:2', '2,4,6'
# - Exercise ambiguous resid handling when chain omitted in AA: if multiple chains share same resid, the first chain (lexicographically) is chosen with an info message.