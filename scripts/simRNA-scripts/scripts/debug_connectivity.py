#!/usr/bin/env python3
"""
Debug script to examine connected and disconnected structures.
Saves examples to files for manual inspection.
"""
import argparse
import os
import sys
from collections import defaultdict
from typing import List, Tuple, Dict

# Import from analyze_strand_connectivity
from analyze_strand_connectivity import (
    parse_ss_detected_file,
    get_strand_for_position,
    extract_base_pairs,
    are_all_strands_connected
)

def analyze_structure_detailed(ss_file: str, temp_level: int = None) -> Dict:
    """
    Analyze a structure in detail, returning all connectivity information.
    
    Returns:
        Dictionary with detailed analysis results
    """
    try:
        structure, strand_ends = parse_ss_detected_file(ss_file)
        num_strands = len(strand_ends) + 1
        pairs = extract_base_pairs(structure)
        
        # Map pairs to strands
        inter_strand_pairs = []
        intra_strand_pairs = []
        strand_connections = defaultdict(set)
        
        for pos1, pos2 in pairs:
            strand1 = get_strand_for_position(pos1, strand_ends)
            strand2 = get_strand_for_position(pos2, strand_ends)
            
            if strand1 != strand2:
                inter_strand_pairs.append((pos1, pos2, strand1, strand2))
                strand_connections[strand1].add(strand2)
                strand_connections[strand2].add(strand1)
            else:
                intra_strand_pairs.append((pos1, pos2, strand1))
        
        is_connected = are_all_strands_connected(pairs, strand_ends)
        
        # Extract frame number from filename
        frame_num = None
        basename = os.path.basename(ss_file)
        import re
        frame_match = re.search(r'-(\d{6})\.ss_detected', basename)
        if frame_match:
            frame_num = int(frame_match.group(1))
        
        return {
            'file': ss_file,
            'frame': frame_num,
            'temperature': temp_level,
            'structure': structure,
            'structure_with_spaces': structure_with_spaces_restored(structure, strand_ends),
            'num_strands': num_strands,
            'strand_ends': strand_ends,
            'total_pairs': len(pairs),
            'intra_strand_pairs': len(intra_strand_pairs),
            'inter_strand_pairs': len(inter_strand_pairs),
            'inter_strand_pairs_list': inter_strand_pairs,
            'strand_connections': dict(strand_connections),
            'is_connected': is_connected
        }
    except Exception as e:
        print(f"Error analyzing {ss_file}: {e}")
        return None

def structure_with_spaces_restored(structure: str, strand_ends: List[int]) -> str:
    """Restore spaces in the structure string at strand boundaries"""
    if not strand_ends:
        return structure
    
    # Handle multiline structures - use only first line for display
    lines = structure.split('\n') if '\n' in structure else [structure]
    first_line = lines[0].strip()
    
    result = []
    prev_end = 0
    for end in strand_ends:
        result.append(first_line[prev_end:end])
        prev_end = end
    result.append(first_line[prev_end:])
    
    return ' '.join(result)

def format_connectivity_report(analysis: Dict) -> str:
    """Format a detailed connectivity report for a structure"""
    lines = []
    lines.append("=" * 80)
    lines.append(f"File: {os.path.basename(analysis['file'])}")
    if analysis['frame']:
        lines.append(f"Frame: {analysis['frame']}")
    if analysis['temperature']:
        lines.append(f"Temperature: {analysis['temperature']}")
    lines.append(f"Status: {'CONNECTED' if analysis['is_connected'] else 'DISCONNECTED'}")
    lines.append("-" * 80)
    
    lines.append(f"Structure: {analysis['structure_with_spaces']}")
    lines.append(f"Number of strands: {analysis['num_strands']}")
    lines.append(f"Strand boundaries (0-based): {analysis['strand_ends']}")
    lines.append("")
    
    lines.append(f"Base pairs:")
    lines.append(f"  Total: {analysis['total_pairs']}")
    lines.append(f"  Intra-strand: {analysis['intra_strand_pairs']}")
    lines.append(f"  Inter-strand: {analysis['inter_strand_pairs']}")
    lines.append("")
    
    if analysis['inter_strand_pairs_list']:
        lines.append("Inter-strand pairs (pos1, pos2, strand1, strand2):")
        for pos1, pos2, s1, s2 in analysis['inter_strand_pairs_list']:
            lines.append(f"  Position {pos1} (strand {s1}) <-> Position {pos2} (strand {s2})")
        lines.append("")
    
    lines.append("Strand connectivity graph:")
    for strand, connections in sorted(analysis['strand_connections'].items()):
        lines.append(f"  Strand {strand} connects to: {sorted(connections)}")
    
    if not analysis['is_connected']:
        lines.append("")
        lines.append("REASON FOR DISCONNECTION:")
        connected_strands = set(analysis['strand_connections'].keys())
        all_strands = set(range(analysis['num_strands']))
        isolated_strands = all_strands - connected_strands
        
        if isolated_strands:
            lines.append(f"  Isolated strands (no inter-strand pairs): {sorted(isolated_strands)}")
        
        # Check for disconnected components
        if len(connected_strands) > 0:
            visited = set()
            start = next(iter(connected_strands))
            to_visit = [start]
            
            while to_visit:
                current = to_visit.pop()
                if current in visited:
                    continue
                visited.add(current)
                
                for neighbor in analysis['strand_connections'][current]:
                    if neighbor not in visited:
                        to_visit.append(neighbor)
            
            if len(visited) < len(connected_strands):
                unreached = connected_strands - visited
                lines.append(f"  Disconnected components: Component 1={sorted(visited)}, Component 2={sorted(unreached)}")
    
    lines.append("=" * 80)
    lines.append("")
    
    return "\n".join(lines)

def collect_examples(temp_level: int, pos_suffix: str = "", 
                     num_connected: int = 5, num_disconnected: int = 5):
    """
    Collect example structures from a temperature level.
    
    Returns:
        tuple: (list of connected examples, list of disconnected examples)
    """
    # Find .ss_detected files for this temperature
    ss_files = []
    
    # Check in pdbs directory
    if os.path.exists("pdbs"):
        ss_files.extend([
            os.path.join("pdbs", f) 
            for f in os.listdir("pdbs") 
            if f.endswith(".ss_detected")
        ])
    
    # Check for level-specific files
    level_pattern = f"level{temp_level}_temp_frames{pos_suffix}"
    for f in os.listdir("."):
        if f.startswith(level_pattern) and f.endswith(".ss_detected"):
            base = f.replace(".ss_detected", "")
            i = 1
            while True:
                ss_file = f"{base}-{i:06d}.ss_detected"
                if os.path.exists(ss_file):
                    ss_files.append(ss_file)
                    i += 1
                else:
                    break
    
    connected_examples = []
    disconnected_examples = []
    
    for ss_file in ss_files:
        if len(connected_examples) >= num_connected and len(disconnected_examples) >= num_disconnected:
            break
        
        analysis = analyze_structure_detailed(ss_file, temp_level)
        if not analysis:
            continue
        
        if analysis['is_connected'] and len(connected_examples) < num_connected:
            connected_examples.append(analysis)
        elif not analysis['is_connected'] and len(disconnected_examples) < num_disconnected:
            disconnected_examples.append(analysis)
    
    return connected_examples, disconnected_examples

def main():
    parser = argparse.ArgumentParser(
        description='Debug strand connectivity analysis'
    )
    parser.add_argument('--name', '-n', required=True,
                       help='Base name for files')
    parser.add_argument('--min-temp', '-t', type=int, required=True,
                       help='Minimum temperature level')
    parser.add_argument('--max-temp', '-T', type=int, required=True,
                       help='Maximum temperature level')
    parser.add_argument('--positions-suffix', default="",
                       help='Position suffix (e.g., "_pos_24to29")')
    parser.add_argument('--num-examples', type=int, default=5,
                       help='Number of examples per category (default: 5)')
    parser.add_argument('--output-dir', default="connectivity_debug",
                       help='Output directory for debug files')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Debugging connectivity analysis: {args.name}")
    print(f"Temperature range: {args.min_temp}-{args.max_temp}")
    print(f"Output directory: {args.output_dir}")
    print("")
    
    # Process each temperature level
    for temp_level in range(args.min_temp, args.max_temp + 1):
        print(f"Processing temperature level {temp_level}...")
        
        connected_examples, disconnected_examples = collect_examples(
            temp_level, 
            args.positions_suffix, 
            args.num_examples, 
            args.num_examples
        )
        
        # Write connected examples
        if connected_examples:
            output_file = os.path.join(
                args.output_dir, 
                f"{args.name}_TEMP{temp_level}_connected{args.positions_suffix}.txt"
            )
            with open(output_file, 'w') as f:
                f.write(f"CONNECTED STRUCTURES - Temperature {temp_level}\n")
                f.write(f"Found {len(connected_examples)} examples\n")
                f.write("\n")
                
                for analysis in connected_examples:
                    f.write(format_connectivity_report(analysis))
            
            print(f"  Saved {len(connected_examples)} connected examples to {output_file}")
        else:
            print(f"  No connected structures found")
        
        # Write disconnected examples
        if disconnected_examples:
            output_file = os.path.join(
                args.output_dir, 
                f"{args.name}_TEMP{temp_level}_disconnected{args.positions_suffix}.txt"
            )
            with open(output_file, 'w') as f:
                f.write(f"DISCONNECTED STRUCTURES - Temperature {temp_level}\n")
                f.write(f"Found {len(disconnected_examples)} examples\n")
                f.write("\n")
                
                for analysis in disconnected_examples:
                    f.write(format_connectivity_report(analysis))
            
            print(f"  Saved {len(disconnected_examples)} disconnected examples to {output_file}")
        else:
            print(f"  No disconnected structures found")
        
        print("")
    
    print("✓ Debug analysis complete")
    print(f"Check the '{args.output_dir}' directory for detailed reports")

if __name__ == "__main__":
    main()