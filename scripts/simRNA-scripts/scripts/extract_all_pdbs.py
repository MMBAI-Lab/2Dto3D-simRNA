#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
from pathlib import Path

def run_simrna_commands(directory, name, suffixes, specific_trafl=None):
    """
    Run SimRNA_trafl2pdbs commands for the specified name, directory, and suffixes.
    """
    # Change to the specified directory
    original_dir = os.getcwd()
    os.chdir(directory)
    
    try:
        # Create pdbs directory if it doesn't exist
        pdbs_dir = Path("pdbs")
        pdbs_dir.mkdir(exist_ok=True)
        
        # Define the reference PDB file
        reference_pdb = f"{name}_01-000001.pdb"
        
        # Check if reference PDB exists
        if not os.path.exists(reference_pdb):
            print(f"Error: Reference PDB file {reference_pdb} not found in {directory}")
            return False
        
        # Define the trajectory files to process
        if specific_trafl:
            # Process only the specific file
            trafl_files = [specific_trafl]
        else:
            # Process files based on suffixes
            trafl_files = [f"{name}_{suffix}.trafl" for suffix in suffixes]
        
        # Check which trajectory files exist
        existing_trafl = []
        for trafl in trafl_files:
            if os.path.exists(trafl):
                existing_trafl.append(trafl)
            else:
                print(f"Warning: {trafl} not found, skipping...")
        
        if not existing_trafl:
            print("Error: No trajectory files found to process")
            return False
        
        print(f"Found reference PDB: {reference_pdb}")
        print(f"Found trajectory files: {existing_trafl}")
        print(f"Output directory: {pdbs_dir.absolute()}")
        print("-" * 60)
        
        # Run SimRNA_trafl2pdbs for each trajectory file
        success_count = 0
        for trafl in existing_trafl:
            print(f"Processing {trafl}...")
            
            # Construct the command
            cmd = ["SimRNA_trafl2pdbs", reference_pdb, trafl, ":"]
            
            try:
                # Run the command and capture output
                result = subprocess.run(
                    cmd, 
                    capture_output=True, 
                    text=True, 
                    check=True,
                    cwd=directory
                )
                
                print(f"✓ Successfully processed {trafl}")
                if result.stdout:
                    print(f"  Output: {result.stdout.strip()}")
                
                # Move generated PDB files to pdbs directory
                # SimRNA typically generates files with pattern: {name}_XX-XXXXXX.pdb
                trafl_base = trafl.replace('.trafl', '')
                pattern = f"{trafl_base}-*.pdb"
                
                # Find and move generated PDB files (excluding starting frames)
                import glob
                generated_pdbs = glob.glob(pattern)
                
                if generated_pdbs:
                    moved_count = 0
                    for pdb_file in generated_pdbs:
                        # Skip starting frames (files ending in 000001.pdb)
                        if pdb_file.endswith('-000001.pdb'):
                            print(f"  Skipping starting frame: {pdb_file}")
                            continue
                        
                        dest_path = pdbs_dir / pdb_file
                        os.rename(pdb_file, dest_path)
                        print(f"  Moved {pdb_file} to {dest_path}")
                        moved_count += 1
                    
                    if moved_count == 0:
                        print(f"  No PDB files to move (only starting frames found)")
                else:
                    print(f"  Warning: No PDB files found with pattern {pattern}")
                
                # Move generated .ss_detected files to pdbs directory (excluding starting frames)
                ss_pattern = f"{trafl_base}-*.ss_detected"
                generated_ss = glob.glob(ss_pattern)
                
                if generated_ss:
                    moved_ss_count = 0
                    for ss_file in generated_ss:
                        # Skip starting frames (files ending in 000001.ss_detected)
                        if ss_file.endswith('-000001.ss_detected'):
                            print(f"  Skipping starting frame: {ss_file}")
                            continue
                        
                        dest_path = pdbs_dir / ss_file
                        os.rename(ss_file, dest_path)
                        print(f"  Moved {ss_file} to {dest_path}")
                        moved_ss_count += 1
                    
                    if moved_ss_count == 0:
                        print(f"  No .ss_detected files to move (only starting frames found)")
                else:
                    print(f"  Warning: No .ss_detected files found with pattern {ss_pattern}")
                
                success_count += 1
                
            except subprocess.CalledProcessError as e:
                print(f"✗ Error processing {trafl}:")
                print(f"  Command failed with return code {e.returncode}")
                if e.stderr:
                    print(f"  Error output: {e.stderr.strip()}")
                continue
            
            except FileNotFoundError:
                print(f"✗ Error: SimRNA_trafl2pdbs command not found")
                print("  Make sure SimRNA is installed and in your PATH")
                return False
        
        print("-" * 60)
        print(f"Processing complete: {success_count}/{len(existing_trafl)} files processed successfully")
        
        # List contents of pdbs directory
        pdb_files = list(pdbs_dir.glob("*.pdb"))
        ss_files = list(pdbs_dir.glob("*.ss_detected"))
        
        if pdb_files:
            print(f"Generated {len(pdb_files)} PDB files in {pdbs_dir}:")
            for pdb in sorted(pdb_files):
                print(f"  {pdb.name}")
        
        if ss_files:
            print(f"Generated {len(ss_files)} .ss_detected files in {pdbs_dir}:")
            for ss in sorted(ss_files):
                print(f"  {ss.name}")
        
        if not pdb_files and not ss_files:
            print("No PDB or .ss_detected files found in output directory")
        
        return success_count > 0
        
    finally:
        # Return to original directory
        os.chdir(original_dir)

def main():
    parser = argparse.ArgumentParser(description='Run SimRNA_trafl2pdbs commands for trajectory files')
    parser.add_argument('directory', help='Directory containing the trajectory files')
    parser.add_argument('--name', '-n', required=True, help='Base name for the files (e.g., "structure")')
    parser.add_argument('--suffixes', '-s', default=','.join([f'{i:02d}' for i in range(1, 21)]),
                        help='Comma-separated list of trajectory suffixes to process (default: 01,02,...,20)')
    parser.add_argument('--trafl', '-t', help='Specific .trafl file to process (overrides --suffixes)')
    parser.add_argument('--dry-run', action='store_true', help='Show commands that would be run without executing them')
    
    args = parser.parse_args()
    
    # Validate directory
    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a valid directory.")
        sys.exit(1)
    
    directory_path = os.path.abspath(args.directory)
    
    # Handle specific trafl file
    specific_trafl = None
    if args.trafl:
        specific_trafl = args.trafl
        # Check if the specific file exists
        trafl_path = os.path.join(directory_path, specific_trafl)
        if not os.path.exists(trafl_path):
            print(f"Error: Specific trafl file {specific_trafl} not found in {directory_path}")
            sys.exit(1)
        print(f"Processing specific file: {specific_trafl}")
    else:
        # Parse suffixes
        suffixes = [s.zfill(2) for s in args.suffixes.split(',') if s.strip()]
        print(f"Processing files with suffixes: {', '.join(suffixes)}")
    
    if args.dry_run:
        print("DRY RUN MODE - Commands that would be executed:")
        print(f"Working directory: {directory_path}")
        print(f"Output directory: {directory_path}/pdbs")
        print()
        
        reference_pdb = f"{args.name}_01-000001.pdb"
        
        if specific_trafl:
            print(f"SimRNA_trafl2pdbs {reference_pdb} {specific_trafl} :")
        else:
            suffixes = [s.zfill(2) for s in args.suffixes.split(',') if s.strip()]
            trafl_files = [f"{args.name}_{suffix}.trafl" for suffix in suffixes]
            for trafl in trafl_files:
                print(f"SimRNA_trafl2pdbs {reference_pdb} {trafl} :")
        
        print("\nOutput files (PDB and .ss_detected) would be moved to pdbs/ directory")
        print("Note: Starting frames (files ending in 000001.pdb and 000001.ss_detected) will be skipped")
        return

    print(f"Running SimRNA trajectory processing...")
    print(f"Directory: {directory_path}")
    print(f"Base name: {args.name}")
    
    if specific_trafl:
        print(f"Specific file: {specific_trafl}")
        suffixes = None
    else:
        suffixes = [s.zfill(2) for s in args.suffixes.split(',') if s.strip()]
        print(f"Suffixes: {', '.join(suffixes)}")
    
    print()
    
    success = run_simrna_commands(directory_path, args.name, suffixes, specific_trafl)
    
    if success:
        print("\n✓ SimRNA processing completed successfully")
    else:
        print("\n✗ SimRNA processing failed")
        sys.exit(1)

if __name__ == "__main__":
    main()