#!/usr/bin/env python3
"""
sasa_from_simrna.py

Compute SASA for specific nucleotide(s) from SimRNA output (.trafl + reference PDB).

Usage (short):
    python sasa_from_simrna.py --pdb reference.pdb --trafl traj.trafl --nucleotide A:10 --aa --probe-radius 1.4 --out per_frame.csv --plot out_box.png
"""

from __future__ import annotations
import argparse
import os
import sys
import tempfile
import uuid
from typing import List, Optional
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Import shared utilities
from _simrna_helpers import (
    CHECK, CROSS, INFO, WARN,
    SimRNAError, SelectorError, TraflParseError,
    detect_simrna_trafl2pdbs, run_simrna_trafl2pdbs,
    parse_selector, parse_frames_arg, format_selector,
    parse_trafl_frames, get_cg_bead_index, validate_cg_trajectory
)


# ----------------------- Shrake-Rupley for beads (approx) -----------------------

def golden_sphere_points(n: int) -> np.ndarray:
    """Return n points approximately evenly distributed on unit sphere (shape (n,3))."""
    indices = np.arange(0, n, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/n)
    theta = math.pi * (1 + 5**0.5) * indices
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    pts = np.vstack([x, y, z]).T
    return pts


def shrake_rupley_beads_frame(coords: np.ndarray, radii: np.ndarray, probe: float = 1.4, n_sphere_points: int = 300, show_progress: bool = False) -> np.ndarray:
    """Compute SASA per bead for a single frame using Shrake-Rupley sampling.
    coords: (N,3), radii: (N,) in Å. probe in Å.
    Returns areas array length N in Å^2.
    """
    N = coords.shape[0]
    pts_unit = golden_sphere_points(n_sphere_points)  # (M,3)
    areas = np.zeros(N, dtype=float)
    radii = radii.astype(float)
    
    iterator = tqdm(range(N), desc="Computing SASA per bead") if (show_progress) else range(N)
    
    for i in iterator:
        Ri = radii[i]
        sphere_rad = Ri + probe
        pts = coords[i] + pts_unit * sphere_rad
        exposed = np.ones(len(pts), dtype=bool)
        for j in range(N):
            if j == i:
                continue
            rj = radii[j]
            d2 = np.sum((pts - coords[j])**2, axis=1)
            exposed &= (d2 > (rj + 1e-8)**2)
            if not np.any(exposed):
                break
        frac = np.count_nonzero(exposed) / float(len(pts))
        areas[i] = frac * (4.0 * math.pi * sphere_rad * sphere_rad)
    return areas


# default bead radii mapping (Å) for the 5 SimRNA beads per residue: P, C4', N1/N9, B1, B2
DEFAULT_BEAD_RADII = {
    'P': 1.8,
    "C4'": 1.7,
    'N': 1.6,
    'B1': 1.7,
    'B2': 1.7,
}


def default_radii_array(n_residues: int) -> np.ndarray:
    """Create default radii array for CG beads."""
    names = ['P', "C4'", 'N', 'B1', 'B2']
    radii = np.array([DEFAULT_BEAD_RADII[n] for n in names] * n_residues, dtype=float)
    return radii


# ----------------------- AA SASA calculator (FreeSASA only) -----------------------

def compute_sasa_freesasa_for_pdb(pdb_path: str, selection: str, probe: float) -> float:
    """
    Compute SASA for `selection` using FreeSASA.
    `selection` should be a selection string like 'chain A and resi 10' or 'resi 10'.
    Returns SASA in Å^2 (sum over atoms in selection).
    """

    structure = freesasa.Structure(pdb_path)
    params = freesasa.Parameters()
    params.setProbeRadius(probe)
    result = freesasa.calc(structure, params)
    selname = 's1'
    selcmd = f"{selname}, {selection}"
    areas = freesasa.selectArea([selcmd], structure, result)
    return float(areas[selname])


def compute_reference_sasa_freesasa(pdb_path: str, selection: str, probe: float) -> float:
    """
    Compute reference SASA for a single residue in isolation.
    This creates a temporary PDB with only the selected residue to get its maximum SASA.
    """
    
    # Read the PDB and extract only the selected residue
    temp_pdb_lines = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                # Parse selection - this is a simplified parser
                if 'chain' in selection and 'resi' in selection:
                    # Extract chain and residue from selection like "chain A and resi 10"
                    parts = selection.split()
                    chain_idx = parts.index('chain') + 1
                    resi_idx = parts.index('resi') + 1
                    target_chain = parts[chain_idx]
                    target_resi = int(parts[resi_idx])
                    
                    line_chain = line[21]
                    line_resi = int(line[22:26].strip())
                    
                    if line_chain == target_chain and line_resi == target_resi:
                        temp_pdb_lines.append(line)
                elif 'resi' in selection:
                    # Extract residue from selection like "resi 10"
                    parts = selection.split()
                    resi_idx = parts.index('resi') + 1
                    target_resi = int(parts[resi_idx])
                    
                    line_resi = int(line[22:26].strip())
                    
                    if line_resi == target_resi:
                        temp_pdb_lines.append(line)
    
    if not temp_pdb_lines:
        return 0.0
    
    # Create temporary PDB file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_f:
        for line in temp_pdb_lines:
            tmp_f.write(line)
        tmp_pdb_path = tmp_f.name
    
    try:
        # Calculate SASA for the isolated residue
        structure = freesasa.Structure(tmp_pdb_path)
        params = freesasa.Parameters()
        params.setProbeRadius(probe)
        result = freesasa.calc(structure, params)
        
        # For isolated residue, just get the total area (all atoms)
        total_sasa = result.totalArea()
        return float(total_sasa)
        
    finally:
        # Clean up temporary file
        try:
            os.unlink(tmp_pdb_path)
        except:
            pass


# ----------------------- High-level driver -----------------------

def run_analysis(args):
    """Main analysis function."""
    # Validate inputs
    if not os.path.isfile(args.trafl):
        raise FileNotFoundError(f"trafl file not found: {args.trafl}")
    if not os.path.isfile(args.pdb):
        raise FileNotFoundError(f"reference PDB not found: {args.pdb}")

    # Parse positions using helper function
    positions = []
    if args.nucleotide:
        tokens = args.nucleotide.split(',')
        for token in tokens:
            try:
                selector = parse_selector(token.strip())
                positions.append(selector)
            except SelectorError as e:
                raise ValueError(f"Invalid nucleotide selector '{token}': {e}")
    else:
        raise ValueError('Please specify --nucleotide (e.g. A:10 or 10)')

    per_frame_records = []
    reference_sasa_cache = {}  # Cache reference SASA values
    
    if args.aa:
        simrna_bin = detect_simrna_trafl2pdbs()
        if not simrna_bin:
            print(f'{INFO} SimRNA_trafl2pdbs not found; falling back to CG mode.')
            args.aa = False
        else:
            try:
                # Parse frames specification
                frames_list = None
                if args.frames and args.frames.strip() != ":":
                    frames_list = parse_frames_arg(args.frames)

                # Run SimRNA conversion
                workdir, pdb_frames = run_simrna_trafl2pdbs(
                    trafl=args.trafl,
                    pdb_first_frame=args.pdb,
                    simrna_bin=simrna_bin,
                    frames=frames_list,
                    aa_mode=True
                )

                if len(pdb_frames) == 0:
                    print(f'{INFO} No PDB frames produced by converter; falling back to CG mode')
                    args.aa = False
                else:
                    print(f"{CHECK} Converted {len(pdb_frames)} frames to all-atom PDBs")
                    
                    # Pre-compute reference SASA values
                    print("Computing reference SASA values for isolated residues...")
                    for selector in positions:
                        chain, resid = selector
                        if chain:
                            selection = f"chain {chain} and resi {resid}"
                        else:
                            selection = f"resi {resid}"
                        key = format_selector(chain, resid)
                        
                        # Use first frame to compute reference SASA
                        first_pdb = sorted(pdb_frames)[0]
                        try:
                            ref_sasa = compute_reference_sasa_freesasa(str(first_pdb), selection, probe=args.probe)
                            reference_sasa_cache[key] = ref_sasa
                            print(f"  {key}: {ref_sasa:.2f} Ų (reference)")
                        except Exception as e:
                            print(f"  Warning: Could not compute reference SASA for {key}: {e}")
                            reference_sasa_cache[key] = None
                    
                    # Process PDB files with progress bar
                    print(f"Computing SASA for {len(positions)} nucleotide(s) across {len(pdb_frames)} frames...")
                    
                    frame_iterator = tqdm(enumerate(sorted(pdb_frames)), total=len(pdb_frames), desc="Processing frames")
                    
                    for fi, pdbf in frame_iterator:
                        for selector in positions:
                            chain, resid = selector
                            if chain:
                                selection = f"chain {chain} and resi {resid}"
                            else:
                                selection = f"resi {resid}"
                            key = format_selector(chain, resid)
                            
                            try:
                                sasa_val = compute_sasa_freesasa_for_pdb(str(pdbf), selection, probe=args.probe)
                                
                                # Calculate relative SASA
                                ref_sasa = reference_sasa_cache.get(key)
                                if ref_sasa and ref_sasa > 0:
                                    rel_sasa = sasa_val / ref_sasa
                                else:
                                    rel_sasa = float('nan')
                                    
                            except Exception as e:
                                sasa_val = float('nan')
                                rel_sasa = float('nan')
                            
                            per_frame_records.append(dict(
                                frame=fi+1, 
                                nucleotide=key, 
                                sasa_A2=sasa_val, 
                                rel_sasa=rel_sasa,
                                method='freesasa', 
                                probe_radius=args.probe
                            ))
                            
                        # Remove pdb if it was created in tempdir and user didn't ask to keep
                        if not args.keep_pdb:
                            try:
                                pdbf.unlink()
                            except Exception:
                                pass
                    
                    print(f"{CHECK} All-atom SASA computation complete")
                    
                    # Clean up temporary directory
                    if not args.keep_pdb and workdir.exists():
                        try:
                            workdir.rmdir()
                        except Exception:
                            pass

            except SimRNAError as e:
                print(f'{INFO} AA conversion failed: {e}. Falling back to CG mode.')
                args.aa = False

    if not args.aa:
        # CG fallback with progress
        print("Parsing trajectory file...")
        headers, coords_array, n_frames, n_atoms = parse_trafl_frames(args.trafl)
        
        # Validate CG trajectory and selectors
        n_residues = validate_cg_trajectory(n_atoms, positions)
        print(f"{CHECK} Found {n_frames} frames, {n_residues} residues (coarse-grained)")
        
        # Handle frame selection
        if args.frames and args.frames.strip() != ":":
            frames_list = parse_frames_arg(args.frames, total_frames_hint=n_frames)
            # Convert to 0-based indices
            frame_indices = [i - 1 for i in frames_list]
            coords_array = coords_array[frame_indices]
            n_frames = coords_array.shape[0]
        
        bead_radii = default_radii_array(n_residues)
        
        print(f"Computing SASA for {len(positions)} nucleotide(s) using coarse-grained beads...")
        print("Note: Relative SASA not available for coarse-grained mode")
        
        frame_iterator = tqdm(range(n_frames), desc="Processing frames")
        
        for fi in frame_iterator:
            frame_coords = coords_array[fi]
            # Show progress for bead calculations only for large systems
            show_bead_progress = (n_atoms > 1000 and fi == 0)  # Only show for first frame of large systems
            areas = shrake_rupley_beads_frame(
                frame_coords, bead_radii, probe=args.probe, 
                n_sphere_points=args.n_sphere_points, 
                show_progress=show_bead_progress
            )
            
            for selector in positions:
                chain, resid = selector
                # Get indices for all 5 beads of this residue
                res_start_idx = (resid - 1) * 5
                bead_slice = slice(res_start_idx, res_start_idx + 5)
                sasa_res = float(np.sum(areas[bead_slice]))
                
                per_frame_records.append(dict(
                    frame=fi+1, 
                    nucleotide=format_selector(None, resid), 
                    sasa_A2=sasa_res, 
                    rel_sasa=float('nan'),  # Not available for CG
                    method='cg_beads', 
                    probe_radius=args.probe
                ))
        
        print(f"{CHECK} Coarse-grained SASA computation complete")

    # Save results
    if len(per_frame_records) > 0:
        df = pd.DataFrame(per_frame_records)
        df.to_csv(args.out, index=False)
        print(f"{CHECK} Saved per-frame SASA to {args.out}")
        
        # plotting + aggregation - now with both absolute and relative SASA
        if args.plot:
            print("Creating boxplots...")
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
            
            # Absolute SASA plot
            groups_abs = df.groupby('nucleotide')['sasa_A2'].apply(list)
            labels = list(groups_abs.index)
            data_abs = [groups_abs[l] for l in labels]
            ax1.boxplot(data_abs)
            ax1.set_xticks(np.arange(1, len(labels)+1))
            ax1.set_xticklabels(labels)
            ax1.set_ylabel('SASA (Ų)')
            ax1.set_title('Absolute SASA distribution per nucleotide')
            ax1.grid(True, alpha=0.3)
            
            # Relative SASA plot
            if not df['rel_sasa'].isna().all():
                groups_rel = df.groupby('nucleotide')['rel_sasa'].apply(list)
                data_rel = [groups_rel[l] for l in labels]
                ax2.boxplot(data_rel)
                ax2.set_xticks(np.arange(1, len(labels)+1))
                ax2.set_xticklabels(labels)
                ax2.set_ylabel('Relative SASA')
                ax2.set_title('Relative SASA distribution per nucleotide')
                ax2.axhline(y=1.0, color='r', linestyle='--', alpha=0.7, label='Fully exposed')
                ax2.legend()
                ax2.grid(True, alpha=0.3)
            else:
                ax2.text(0.5, 0.5, 'Relative SASA not available\n(coarse-grained mode)', 
                        ha='center', va='center', transform=ax2.transAxes)
                ax2.set_title('Relative SASA distribution per nucleotide')
            
            plt.tight_layout()
            plt.savefig(args.plot, dpi=300, bbox_inches='tight')
            print(f"{CHECK} Saved boxplots to {args.plot}")
            
        if args.append_run_rows:
            run_id = args.run_id or str(uuid.uuid4())
            df2 = df.copy()
            df2['run_id'] = run_id
            df2['method'] = df2.get('method', 'unknown')
            df2['probe_radius'] = args.probe
            df2 = df2[['run_id','nucleotide','method','probe_radius','frame','sasa_A2','rel_sasa']]
            mode = 'a' if os.path.exists(args.append_run_rows) else 'w'
            header = not os.path.exists(args.append_run_rows)
            df2.to_csv(args.append_run_rows, mode=mode, header=header, index=False)
            print(f"{CHECK} Appended results to aggregate CSV: {args.append_run_rows}")
            
        # Update summary to include relative SASA
        summary_abs = df.groupby('nucleotide')['sasa_A2'].agg(['median','mean','std','min','max', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)])
        summary_abs.columns = ['median_abs','mean_abs','std_abs','min_abs','max_abs','q1_abs','q3_abs']
        
        if not df['rel_sasa'].isna().all():
            summary_rel = df.groupby('nucleotide')['rel_sasa'].agg(['median','mean','std','min','max', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)])
            summary_rel.columns = ['median_rel','mean_rel','std_rel','min_rel','max_rel','q1_rel','q3_rel']
            summary = pd.concat([summary_abs, summary_rel], axis=1)
        else:
            summary = summary_abs
        
        summary = summary.reset_index()
        summary['run_id'] = args.run_id or str(uuid.uuid4())
        
        # Write boxplot summary with append logic
        boxplot_csv_path = args.boxplot_csv or (os.path.splitext(args.out)[0]+'.boxplot_summary.csv')
        if os.path.exists(boxplot_csv_path):
            # Append to existing file
            summary.to_csv(boxplot_csv_path, mode='a', header=False, index=False)
            print(f'{CHECK} Appended boxplot summary to existing file: {boxplot_csv_path}')
        else:
            # Create new file with header
            summary.to_csv(boxplot_csv_path, index=False)
            print(f'{CHECK} Created new boxplot summary CSV: {boxplot_csv_path}')


# ----------------------- CLI -----------------------

def build_parser():
    """Build command-line argument parser."""
    p = argparse.ArgumentParser(description='Compute SASA for nucleotides from SimRNA trafl + reference PDB')
    
    # Required arguments
    p.add_argument('--pdb', required=True, help='Reference PDB (used for AA conversion or mapping)')
    p.add_argument('--trafl', required=True, help='SimRNA .trafl trajectory file')
    p.add_argument('--nucleotide', help='nucleotide(s) to analyze, e.g. A:10 or A:10,B:12 or 10 or 10,12 (1-based indices)')
    
    # Analysis options
    p.add_argument('--probe-radius', type=float, default=1.4, dest='probe', 
                   help='Probe radius in Å for SASA (default 1.4 Å)')
    p.add_argument('--aa', action='store_true', 
                   help='Try all-atom reconstruction via SimRNA_trafl2pdbs')
    p.add_argument('--simrna-bin', help='Path to SimRNA_trafl2pdbs executable')
    p.add_argument('--frames', help='Frames spec for SimRNA_trafl2pdbs (e.g. "1:10" or ":" for all)')
    
    # Output options
    p.add_argument('--out', default='sasa_per_frame.csv', help='Output per-frame CSV path')
    p.add_argument('--plot', help='Save boxplot PNG to this path')
    p.add_argument('--append-run-rows', help='Append per-frame rows into this aggregate CSV')
    p.add_argument('--boxplot-csv', help='Path to write boxplot-summary CSV')
    p.add_argument('--run-id', help='Optional run identifier to store in aggregate CSV')
    p.add_argument('--keep-pdb', action='store_true', 
                   help='Keep PDB files generated by SimRNA_trafl2pdbs (default: remove them)')
    p.add_argument('--n-sphere-points', type=int, default=300, 
                   help='Number of sphere sample points for CG Shrake-Rupley. Default 300.')
    
    return p


def main():
    """Main entry point."""
    parser = build_parser()
    args = parser.parse_args()
    try:
        run_analysis(args)
    except (SimRNAError, SelectorError, TraflParseError, ValueError, FileNotFoundError) as e:
        print(f'{CROSS} Error: {e}')
        sys.exit(1)
    except Exception as e:
        print(f'{CROSS} Unexpected error: {e}')
        sys.exit(1)


if __name__ == '__main__':
    main()
