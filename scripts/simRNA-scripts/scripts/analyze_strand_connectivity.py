#!/usr/bin/env python3
"""
Analyze strand connectivity in RNA structures.
Determines what proportion of structures have all sequences connected
(i.e., the complex includes all strands).
"""
import argparse
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
from typing import List, Tuple, Dict, Set

def parse_ss_detected_file(filepath: str) -> Tuple[str, List[int]]:
    """
    Parse .ss_detected file to extract structure and strand boundaries.
    
    Note: Multi-line structures are collapsed to single line. All lines represent
    the same positions, so we only need to look at the first line for strand boundaries.
    
    Returns:
        tuple: (structure_string, list of strand end positions)
    """
    with open(filepath, 'r') as f:
        content = f.read().strip()
    
    # Split into lines - each line is an alternative representation of the same structure
    lines = content.split('\n')
    
    # Use only the first line to determine strand boundaries
    # (all lines represent the same positions and strands)
    first_line = lines[0].strip()
    
    # Find spaces in the first line (strand separators)
    strand_ends = []
    cumulative_pos = 0
    
    # Split by spaces and track positions
    parts = first_line.split(' ')
    for i, part in enumerate(parts[:-1]):  # Exclude last part
        cumulative_pos += len(part)
        strand_ends.append(cumulative_pos)
        cumulative_pos += 1  # Account for the space itself
    
    # Combine all lines and remove spaces to get pure structure
    # Each character position appears once per line
    structure = content.replace(' ', '').replace('\n', '')
    
    return structure, strand_ends

def get_strand_for_position(pos: int, strand_ends: List[int]) -> int:
    """
    Determine which strand a position belongs to.
    
    Args:
        pos: 0-based position in the structure
        strand_ends: List of positions where each strand ends (exclusive)
    
    Returns:
        Strand number (0-based)
    """
    for strand_num, end_pos in enumerate(strand_ends):
        if pos < end_pos:
            return strand_num
    return len(strand_ends)  # Last strand

def extract_base_pairs(structure: str) -> List[Tuple[int, int]]:
    """
    Extract base pairs from dot-bracket notation.
    Handles multi-line structures by processing each line separately
    and combining the results (each line shows different pairs for the same positions).
    
    Returns:
        List of (pos1, pos2) tuples representing base pairs
    """
    pairs = []
    
    # Structure might be multi-line - each line represents the same positions
    # but shows different base pairs (when structure is too complex for single line)
    lines = structure.split('\n') if '\n' in structure else [structure]
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        stack = []
        for i, char in enumerate(line):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if stack:
                    j = stack.pop()
                    pairs.append((j, i))
    
    # Remove duplicates (shouldn't happen but just in case)
    pairs = list(set(pairs))
    
    return pairs

def are_all_strands_connected(pairs: List[Tuple[int, int]], strand_ends: List[int]) -> bool:
    """
    Determine if all strands are connected through base pairs.
    Uses graph connectivity algorithm.
    
    Args:
        pairs: List of base pair tuples
        strand_ends: List of strand boundary positions
    
    Returns:
        True if all strands are connected, False otherwise
    """
    num_strands = len(strand_ends) + 1
    
    if num_strands <= 1:
        return True  # Single strand is trivially connected
    
    # Build adjacency list of which strands are connected
    strand_connections = defaultdict(set)
    
    for pos1, pos2 in pairs:
        strand1 = get_strand_for_position(pos1, strand_ends)
        strand2 = get_strand_for_position(pos2, strand_ends)
        
        if strand1 != strand2:
            strand_connections[strand1].add(strand2)
            strand_connections[strand2].add(strand1)
    
    # Check connectivity using BFS/DFS
    if not strand_connections:
        return False  # No inter-strand connections
    
    visited = set()
    to_visit = [0]  # Start from first strand
    
    while to_visit:
        current = to_visit.pop()
        if current in visited:
            continue
        visited.add(current)
        
        for neighbor in strand_connections[current]:
            if neighbor not in visited:
                to_visit.append(neighbor)
    
    return len(visited) == num_strands

def analyze_ss_file(ss_file: str) -> Tuple[bool, int, int]:
    """
    Analyze a single .ss_detected file.
    
    Returns:
        tuple: (is_connected, num_strands, num_inter_strand_pairs)
    """
    try:
        structure, strand_ends = parse_ss_detected_file(ss_file)
        num_strands = len(strand_ends) + 1
        pairs = extract_base_pairs(structure)
        
        # Count inter-strand pairs
        inter_strand_pairs = 0
        for pos1, pos2 in pairs:
            strand1 = get_strand_for_position(pos1, strand_ends)
            strand2 = get_strand_for_position(pos2, strand_ends)
            if strand1 != strand2:
                inter_strand_pairs += 1
        
        is_connected = are_all_strands_connected(pairs, strand_ends)
        
        return is_connected, num_strands, inter_strand_pairs
    except Exception as e:
        print(f"Error analyzing {ss_file}: {e}")
        return False, 0, 0

def process_temperature_level(temp_level: int, pos_suffix: str = "") -> Dict:
    """
    Process all .ss_detected files for a temperature level.
    
    Returns:
        Dictionary with analysis results
    """
    # Look for .ss_detected files matching this temperature
    pattern = f"*TEMP{temp_level}_*{pos_suffix}.png"
    
    # Find the corresponding .ss_detected files
    ss_files = []
    
    # Check in pdbs directory (if it exists)
    if os.path.exists("pdbs"):
        ss_files.extend([f for f in os.listdir("pdbs") if f.endswith(".ss_detected")])
        ss_files = [os.path.join("pdbs", f) for f in ss_files]
    
    # Also check in current directory for level-specific files
    level_pattern = f"level{temp_level}_temp_frames{pos_suffix}"
    for f in os.listdir("."):
        if f.startswith(level_pattern) and f.endswith(".ss_detected"):
            # This is a trafl metadata file, get associated ss files
            base = f.replace(".ss_detected", "")
            # Look for numbered ss files
            i = 1
            while True:
                ss_file = f"{base}-{i:06d}.ss_detected"
                if os.path.exists(ss_file):
                    ss_files.append(ss_file)
                    i += 1
                else:
                    break
    
    if not ss_files:
        # Try to find from PNG files what ss_detected files were analyzed
        csv_file = f"accuracy_TEMP{temp_level}_pairs{pos_suffix}.csv"
        if os.path.exists(csv_file):
            df = pd.read_csv(csv_file)
            if 'structure' in df.columns:
                for struct_path in df['structure'].unique():
                    ss_file = struct_path.replace('.pdb', '.ss_detected')
                    if os.path.exists(ss_file):
                        ss_files.append(ss_file)
    
    results = {
        'temperature': temp_level,
        'total_structures': 0,
        'connected_structures': 0,
        'disconnected_structures': 0,
        'connectivity_rate': 0.0,
        'num_strands': 0,
        'avg_inter_strand_pairs': 0.0
    }
    
    if not ss_files:
        print(f"Warning: No .ss_detected files found for temperature {temp_level}")
        return results
    
    connected_count = 0
    total_count = 0
    inter_strand_pairs_list = []
    
    for ss_file in ss_files:
        is_connected, num_strands, inter_pairs = analyze_ss_file(ss_file)
        total_count += 1
        if is_connected:
            connected_count += 1
        inter_strand_pairs_list.append(inter_pairs)
        
        if results['num_strands'] == 0:
            results['num_strands'] = num_strands
    
    results['total_structures'] = total_count
    results['connected_structures'] = connected_count
    results['disconnected_structures'] = total_count - connected_count
    results['connectivity_rate'] = (connected_count / total_count * 100) if total_count > 0 else 0.0
    results['avg_inter_strand_pairs'] = np.mean(inter_strand_pairs_list) if inter_strand_pairs_list else 0.0
    
    return results

def create_connectivity_plots(connectivity_data: List[Dict], name: str, pos_suffix: str = ""):
    """
    Create visualization plots for strand connectivity analysis.
    """
    if not connectivity_data:
        print("No connectivity data to plot")
        return []
    
    df = pd.DataFrame(connectivity_data)
    output_files = []
    
    # Plot 1: Connectivity rate by temperature
    plt.figure(figsize=(10, 6))
    plt.plot(df['temperature'], df['connectivity_rate'], marker='o', linewidth=2, markersize=8)
    plt.xlabel('Temperature Level', fontsize=12)
    plt.ylabel('Connectivity Rate (%)', fontsize=12)
    plt.title(f'Strand Connectivity Rate by Temperature\n({name})', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.ylim(-5, 105)
    
    output_file = f"{name}_connectivity_rate{pos_suffix}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    output_files.append(output_file)
    
    # Plot 2: Connected vs Disconnected structures
    plt.figure(figsize=(12, 6))
    x = np.arange(len(df['temperature']))
    width = 0.35
    
    plt.bar(x - width/2, df['connected_structures'], width, label='Connected', color='green', alpha=0.7)
    plt.bar(x + width/2, df['disconnected_structures'], width, label='Disconnected', color='red', alpha=0.7)
    
    plt.xlabel('Temperature Level', fontsize=12)
    plt.ylabel('Number of Structures', fontsize=12)
    plt.title(f'Connected vs Disconnected Structures\n({name})', fontsize=14, fontweight='bold')
    plt.xticks(x, df['temperature'])
    plt.legend()
    plt.grid(True, alpha=0.3, axis='y')
    
    output_file = f"{name}_connectivity_counts{pos_suffix}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    output_files.append(output_file)
    
    # Plot 3: Average inter-strand pairs
    plt.figure(figsize=(10, 6))
    plt.plot(df['temperature'], df['avg_inter_strand_pairs'], marker='s', linewidth=2, 
             markersize=8, color='purple')
    plt.xlabel('Temperature Level', fontsize=12)
    plt.ylabel('Average Inter-Strand Base Pairs', fontsize=12)
    plt.title(f'Inter-Strand Base Pairs by Temperature\n({name})', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    output_file = f"{name}_inter_strand_pairs{pos_suffix}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    output_files.append(output_file)
    
    return output_files

def save_connectivity_summary(connectivity_data: List[Dict], name: str, pos_suffix: str = ""):
    """
    Save connectivity analysis to CSV file.
    """
    if not connectivity_data:
        return None
    
    df = pd.DataFrame(connectivity_data)
    output_file = f"{name}_connectivity_summary{pos_suffix}.csv"
    df.to_csv(output_file, index=False, float_format='%.2f')
    print(f"Connectivity summary saved to: {output_file}")
    return output_file

def main():
    parser = argparse.ArgumentParser(
        description='Analyze strand connectivity in RNA structures'
    )
    parser.add_argument('--name', '-n', required=True, 
                       help='Base name for files')
    parser.add_argument('--min-temp', '-t', type=int, required=True, 
                       help='Minimum temperature level')
    parser.add_argument('--max-temp', '-T', type=int, required=True, 
                       help='Maximum temperature level')
    parser.add_argument('--positions-suffix', default="", 
                       help='Position suffix (e.g., "_pos_24to29")')
    parser.add_argument('--skip-plots', action='store_true', 
                       help='Skip generating plots')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug mode (save examples of connected/disconnected)')
    parser.add_argument('--debug-examples', type=int, default=5,
                       help='Number of debug examples per category (default: 5)')
    parser.add_argument('--quiet', '-q', action='store_true', 
                       help='Reduce output verbosity')
    
    args = parser.parse_args()
    
    if not args.quiet:
        print(f"Analyzing strand connectivity: {args.name}")
        print(f"Temperature range: {args.min_temp}-{args.max_temp}")
    
    # Process each temperature level
    connectivity_data = []
    
    for temp_level in range(args.min_temp, args.max_temp + 1):
        if not args.quiet:
            print(f"Processing temperature level {temp_level}...")
        
        results = process_temperature_level(temp_level, args.positions_suffix)
        
        if results['total_structures'] > 0:
            connectivity_data.append(results)
            if not args.quiet:
                print(f"  Temperature {temp_level}: {results['connectivity_rate']:.1f}% connected "
                      f"({results['connected_structures']}/{results['total_structures']})")
        else:
            if not args.quiet:
                print(f"  Temperature {temp_level}: No structures found")
    
    if not connectivity_data:
        print("No connectivity data found!")
        return
    
    # Save summary
    save_connectivity_summary(connectivity_data, args.name, args.positions_suffix)
    
    # Create plots
    if not args.skip_plots:
        plot_files = create_connectivity_plots(connectivity_data, args.name, args.positions_suffix)
        if not args.quiet:
            print(f"Created {len(plot_files)} connectivity plots")
    
    # Run debug analysis if requested
    if args.debug:
        if not args.quiet:
            print("\nRunning debug analysis...")
        
        debug_script = os.path.join(os.path.dirname(__file__), "debug_connectivity.py")
        if os.path.exists(debug_script):
            import subprocess
            cmd = [
                "python", debug_script,
                "--name", args.name,
                "--min-temp", str(args.min_temp),
                "--max-temp", str(args.max_temp),
                "--positions-suffix", args.positions_suffix,
                "--num-examples", str(args.debug_examples)
            ]
            subprocess.run(cmd)
        else:
            print(f"Warning: debug_connectivity.py not found at {debug_script}")
    
    # Print summary
    if not args.quiet:
        print("\nConnectivity Summary:")
        print("-" * 60)
        for data in connectivity_data:
            print(f"Temp {data['temperature']}: {data['connectivity_rate']:6.1f}% "
                  f"({data['connected_structures']:4d}/{data['total_structures']:4d}) "
                  f"Avg inter-strand pairs: {data['avg_inter_strand_pairs']:.1f}")

if __name__ == "__main__":
    main()