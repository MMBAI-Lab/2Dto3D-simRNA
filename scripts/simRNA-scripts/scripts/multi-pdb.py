#!/usr/bin/env python3

import glob
import os
import argparse

def create_trajectory_pdb(pdb_files, output_file):
    """
    Create a multi-model PDB file that ChimeraX can read as a trajectory
    """
    with open(output_file, 'w') as outfile:
        model_num = 1
        for pdb_file in sorted(pdb_files):
            with open(pdb_file, 'r') as infile:
                # Write MODEL record
                outfile.write(f"MODEL     {model_num:4d}\n")
                # Copy all ATOM/HETATM records
                for line in infile:
                    if line.startswith(('ATOM', 'HETATM')):
                        outfile.write(line)
                # Write ENDMDL record
                outfile.write("ENDMDL\n")
                model_num += 1
        # Write END record
        outfile.write("END\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine all PDB files in a folder into a multi-model trajectory PDB.")
    parser.add_argument('--folder', '-f', default='.', help='Folder containing PDB files (default: current directory)')
    parser.add_argument('--output', '-o', default='trajectory.pdb', help='Output multi-model PDB file (default: trajectory.pdb)')
    args = parser.parse_args()

    pdb_path = os.path.abspath(args.folder)
    search_pattern = os.path.join(pdb_path, "*.pdb")
    pdb_files = glob.glob(search_pattern)

    print(f"Searching for PDB files in: {pdb_path}")
    print(f"Found {len(pdb_files)} PDB files: {pdb_files}")

    if pdb_files:
        create_trajectory_pdb(pdb_files, args.output)
        print(f"Created {args.output}")
    else:
        print("No PDB files found in the folder!")
