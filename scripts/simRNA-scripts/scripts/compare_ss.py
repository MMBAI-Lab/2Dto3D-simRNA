#!/usr/bin/env python3
import os
import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import re
from collections import defaultdict
import pandas as pd
import csv

def read_structure_file(file_path):
    """Read structure file and extract base pairs from all lines"""
    all_pairs = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('//') and any(char in line for char in '().'):
                clean_line = line.replace(' ', '')
                line_pairs = parse_dot_bracket(clean_line)
                all_pairs.update(line_pairs)
    return all_pairs

def parse_dot_bracket(structure):
    """Parse dot-bracket notation and return base pair mappings"""
    stack = []
    pairs = {}
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i + 1)
        elif char == ')' and stack:
            pairs[stack.pop()] = i + 1
    return pairs

def get_structure_length(file_path):
    """Get total structure length"""
    max_length = 0
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('//') and any(char in line for char in '().'):
                max_length = max(max_length, len(line.replace(' ', '')))
    return max_length

def get_unpaired_positions(pairs, structure_length):
    """Get positions not involved in base pairs"""
    paired_positions = set()
    for open_pos, close_pos in pairs.items():
        paired_positions.update([open_pos, close_pos])
    return set(range(1, structure_length + 1)) - paired_positions

def calculate_structure_score(reference_pairs, predicted_pairs, reference_length, predicted_length, include_dots=False, dots_only=False):
    """Calculate score comparing reference and predicted structures"""
    if not reference_pairs and not include_dots and not dots_only:
        return 0, 0, 0, 0, 0, 0
    
    # Base pair score
    pair_score = 0
    total_pairs = len(reference_pairs)
    if not dots_only:
        for ref_open, ref_close in reference_pairs.items():
            if ref_open in predicted_pairs and predicted_pairs[ref_open] == ref_close:
                pair_score += 1
    
    # Dot score
    dot_score = 0
    total_dots = 0
    if include_dots or dots_only:
        reference_dots = get_unpaired_positions(reference_pairs, reference_length)
        predicted_dots = get_unpaired_positions(predicted_pairs, predicted_length)
        total_dots = len(reference_dots)
        
        for dot_pos in reference_dots:
            if dot_pos <= predicted_length and dot_pos in predicted_dots:
                dot_score += 1.0 if dots_only else 0.5
    
    # Total calculations
    if dots_only:
        total_score = dot_score
        total_possible_score = total_dots * 1.0
    else:
        total_score = pair_score + dot_score
        total_possible_score = total_pairs + (total_dots * 0.5 if include_dots else 0)
    
    return total_score, total_possible_score, pair_score, dot_score, total_pairs, total_dots

def calculate_position_specific_score(reference_pairs, predicted_pairs, reference_length, predicted_length, positions, include_dots=False, dots_only=False):
    """Calculate score for specific positions only"""
    if not positions:
        return 0, 0
    
    correct = 0
    total = len(positions)
    
    # Convert positions to 0-based for internal calculations
    positions_0based = [p - 1 for p in positions]
    
    # Create mapping of position to partner for reference
    ref_partners = {}
    for open_pos, close_pos in reference_pairs.items():
        ref_partners[open_pos - 1] = close_pos - 1  # Convert to 0-based
        ref_partners[close_pos - 1] = open_pos - 1
    
    # Create mapping of position to partner for predicted
    pred_partners = {}
    for open_pos, close_pos in predicted_pairs.items():
        pred_partners[open_pos - 1] = close_pos - 1  # Convert to 0-based
        pred_partners[close_pos - 1] = open_pos - 1
    
    # Check only the specified positions
    for pos in positions_0based:
        if pos < reference_length and pos < predicted_length:
            ref_partner = ref_partners.get(pos, None)
            pred_partner = pred_partners.get(pos, None)
            
            # Check if the pairing is correct
            if ref_partner == pred_partner:
                correct += 1
    
    return correct, total

def extract_trajectory_info(filename):
    """Extract trajectory suffix and frame number from filename"""
    patterns = [
        r'.*_(\d{2})-(\d{6})\.ss_detected$',
        r'all_thrs[\d.]+A_clust(\d+)-(\d{6})\.ss_detected$'
    ]
    
    for pattern in patterns:
        match = re.match(pattern, filename)
        if match:
            return match.group(1).zfill(2), int(match.group(2))
    
    if filename.endswith('.ss_detected'):
        numbers = re.findall(r'\d+', filename)
        if numbers:
            return "01", int(numbers[-1])
    
    return None, None

def process_structures(str_file, ss_folder, include_dots=False, dots_only=False, positions=None, position_threshold=None):
    """Process all .ss_detected files and calculate scores"""
    reference_pairs = read_structure_file(str_file)
    reference_length = get_structure_length(str_file)
    
    print(f"Reference: {reference_length} nucleotides, {len(reference_pairs)} pairs")
    if include_dots or dots_only:
        reference_dots = get_unpaired_positions(reference_pairs, reference_length)
        print(f"Reference has {len(reference_dots)} unpaired positions")
    
    if positions:
        print(f"Position-specific analysis for positions: {positions}")
    
    # Group files by trajectory
    trajectory_groups = {}
    for file in os.listdir(ss_folder):
        if file.endswith('.ss_detected'):
            suffix, frame_num = extract_trajectory_info(file)
            if suffix is not None:
                if suffix not in trajectory_groups:
                    trajectory_groups[suffix] = []
                trajectory_groups[suffix].append({'filename': file, 'frame_number': frame_num})
    
    if not trajectory_groups:
        print(f"No .ss_detected files found in {ss_folder}")
        return {}
    
    # Process each trajectory
    all_results = {}
    total_processed = 0
    total_filtered = 0
    
    for suffix in sorted(trajectory_groups.keys()):
        print(f"\nProcessing trajectory _{suffix}...")
        files_in_group = sorted(trajectory_groups[suffix], key=lambda x: x['frame_number'])
        
        results = []
        for file_info in files_in_group:
            filename = file_info['filename']
            file_path = os.path.join(ss_folder, filename)
            
            try:
                predicted_pairs = read_structure_file(file_path)
                predicted_length = get_structure_length(file_path)
                
                total_score, total_possible, pair_score, dot_score, total_pairs, total_dots = calculate_structure_score(
                    reference_pairs, predicted_pairs, reference_length, predicted_length, include_dots, dots_only
                )
                
                accuracy = total_score / total_possible if total_possible > 0 else 0
                pair_accuracy = pair_score / total_pairs if total_pairs > 0 else 0
                
                if dots_only:
                    dot_accuracy = (dot_score / total_dots) if total_dots > 0 else 0
                else:
                    dot_accuracy = (dot_score / (total_dots * 0.5)) if total_dots > 0 else 0
                
                result = {
                    'file': filename,
                    'frame_number': file_info['frame_number'],
                    'total_score': total_score,
                    'total_possible': total_possible,
                    'pair_score': pair_score,
                    'dot_score': dot_score,
                    'total_pairs': total_pairs,
                    'total_dots': total_dots,
                    'accuracy': accuracy,
                    'pair_accuracy': pair_accuracy,
                    'dot_accuracy': dot_accuracy,
                    'trajectory': suffix
                }
                
                # Calculate position-specific accuracy if positions are specified
                if positions:
                    position_correct, position_total = calculate_position_specific_score(
                        reference_pairs, predicted_pairs, reference_length, predicted_length, positions, include_dots, dots_only
                    )
                    position_accuracy = position_correct / position_total if position_total > 0 else 0
                    result['position_accuracy'] = position_accuracy
                    result['position_correct'] = position_correct
                    result['position_total'] = position_total
                    
                    # Apply position threshold filter if specified
                    if position_threshold is not None and position_accuracy < position_threshold:
                        total_filtered += 1
                        continue  # Skip this frame
                
                if dots_only:
                    print(f"  {filename}: {dot_score:.1f}/{total_dots:.1f} dots (accuracy: {accuracy:.3f})")
                elif include_dots:
                    print(f"  {filename}: {total_score:.1f}/{total_possible:.1f} total (accuracy: {accuracy:.3f})")
                else:
                    print(f"  {filename}: {pair_score}/{total_pairs} pairs (accuracy: {accuracy:.3f})")
                
                if positions:
                    print(f"    Position accuracy: {result['position_accuracy']:.3f} ({result['position_correct']}/{result['position_total']})")
                
                results.append(result)
                total_processed += 1
                
            except Exception as e:
                print(f"  Error processing {filename}: {e}")
        
        all_results[suffix] = results
    
    if position_threshold is not None:
        print(f"\nFiltered {total_filtered} frames below position threshold {position_threshold}")
    print(f"Total processed: {total_processed} frames")
    
    return all_results

def plot_all_results_combined(all_results, output_prefix="rna_structure_comparison", include_dots=False, dots_only=False, positions_only=False):
    """Create single plot with all results ordered by accuracy"""
    if not all_results:
        return
    
    # Combine all results
    all_combined = []
    trajectory_colors = {'02': 'steelblue', '03': 'forestgreen', '04': 'crimson', 
                        '05': 'orange', '06': 'purple', '07': 'brown', '08': 'pink'}
    
    for suffix, results in all_results.items():
        for result in results:
            result_copy = result.copy()
            result_copy['trajectory'] = suffix
            all_combined.append(result_copy)
    
    if not all_combined:
        return
    
    # Choose the appropriate accuracy column
    accuracy_key = 'position_accuracy' if positions_only and 'position_accuracy' in all_combined[0] else 'accuracy'
    
    all_combined.sort(key=lambda x: x[accuracy_key], reverse=True)
    
    # Extract data for plotting
    accuracies = [r[accuracy_key] for r in all_combined]
    trajectories = [r['trajectory'] for r in all_combined]
    
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    
    # Group by trajectory for vectorized plotting
    trajectory_data = {}
    for i, (acc, traj) in enumerate(zip(accuracies, trajectories)):
        if traj not in trajectory_data:
            trajectory_data[traj] = {'x': [], 'y': []}
        trajectory_data[traj]['x'].append(i)
        trajectory_data[traj]['y'].append(acc)
    
    # Plot points
    for traj in sorted(trajectory_data.keys()):
        color = trajectory_colors.get(traj, 'gray')
        ax.scatter(trajectory_data[traj]['x'], trajectory_data[traj]['y'], 
                  c=color, s=50, alpha=0.7, label=f'Trajectory _{traj}')
    
    # Add reference regions
    ax.axhspan(0.9, 1.0, alpha=0.2, color='green', label='90-100% (Excellent)')
    ax.axhspan(0.8, 0.9, alpha=0.2, color='lightgreen', label='80-90% (Good)')
    ax.axhspan(0.7, 0.8, alpha=0.2, color='yellow', label='70-80% (Fair)')
    ax.axhspan(0.6, 0.7, alpha=0.2, color='orange', label='60-70% (Poor)')
    ax.axhspan(0.0, 0.6, alpha=0.2, color='red', label='0-60% (Very Poor)')
    
    # Statistics
    mean_acc = np.mean(accuracies)
    median_acc = np.median(accuracies)
    std_acc = np.std(accuracies)
    
    scoring_method = "dots only" if dots_only else ("with dots" if include_dots else "pairs only")
    stats_text = f'Scoring: {scoring_method}\nMean: {mean_acc:.3f}\nMedian: {median_acc:.3f}\nStd Dev: {std_acc:.3f}\nTotal: {len(all_combined)}'
    ax.text(0.02, 0.02, stats_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Count distribution
    ranges = [(0.9, 1.001, 'Excellent'), (0.8, 0.9, 'Good'), (0.7, 0.8, 'Fair'), 
              (0.6, 0.7, 'Poor'), (0.0, 0.6, 'Very Poor')]
    
    counts_text = "Distribution:\n"
    for low, high, label in ranges:
        count = sum(1 for acc in accuracies if low <= acc < high)
        percentage = (count / len(accuracies)) * 100
        counts_text += f"{label}: {count} ({percentage:.1f}%)\n"
    
    ax.text(0.98, 0.98, counts_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Styling
    position_suffix = " (Position-Specific)" if positions_only else ""
    scoring_method = "dots only" if dots_only else ("with dots" if include_dots else "pairs only")
    
    title_suffix = f" (Dots Only){position_suffix}" if dots_only else (f" (Including Dots){position_suffix}" if include_dots else f" (Base Pairs Only){position_suffix}")
    ax.set_title(f'RNA Structure Comparison: All Results Ordered by Accuracy{title_suffix}', 
                fontsize=20, fontweight='bold', pad=20)
    ax.set_ylabel('Accuracy (Total Score / Total Possible)', fontsize=16, fontweight='bold')
    ax.set_xlabel('Structure Index (Ordered by Accuracy)', fontsize=16, fontweight='bold')
    
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.1%}'))
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    if len(all_combined) > 50:
        step = len(all_combined) // 20
        ax.set_xticks(range(0, len(all_combined), step))
    
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=11, loc='center left', bbox_to_anchor=(1, 0.5))
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    plt.tight_layout()
    
    dots_suffix = "_dots_only" if dots_only else ("_with_dots" if include_dots else "")
    output_file = f"{output_prefix}_all_combined{dots_suffix}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Combined plot saved to {output_file}")
    plt.close()
    
    return output_file

def main():
    parser = argparse.ArgumentParser(description='Compare RNA secondary structures and score base pair accuracy')
    parser.add_argument('str_file', help='Reference .str file containing the target structure')
    parser.add_argument('ss_folder', help='Folder containing .ss_detected files')
    parser.add_argument('--output', '-o', help='Output file prefix for the plots', default='rna_structure_comparison')
    parser.add_argument('--csv', help='Save results to CSV file')
    parser.add_argument('--all', '-all', action='store_true', help='Create combined plot ordered by accuracy')
    parser.add_argument('--dots', action='store_true', help='Include unpaired positions with half weight')
    parser.add_argument('--dots-only', action='store_true', help='Score ONLY unpaired positions')
    parser.add_argument('--positions', help='Comma-separated list of positions to analyze (1-based)')
    parser.add_argument('--position-threshold', type=float, 
                       help='Minimum position accuracy threshold to include frames')
    
    args = parser.parse_args()
    
    if not os.path.isfile(args.str_file) or not os.path.isdir(args.ss_folder):
        print("Error: Invalid file or directory")
        sys.exit(1)
    
    print(f"Comparing against: {args.str_file}")
    print(f"Processing files in: {args.ss_folder}")
    scoring_mode = "dots only" if args.dots_only else ("with dots" if args.dots else "pairs only")
    print(f"Scoring: {scoring_mode}")
    
    # Parse positions if provided
    positions = None
    positions_only = False
    if args.positions:
        try:
            positions = [int(x.strip()) for x in args.positions.split(',')]  # Keep 1-based
            positions_only = True
            print(f"Position-specific analysis for positions: {positions}")
        except ValueError:
            print("Error: Invalid positions format. Use comma-separated numbers like '24,25,26,27'")
            sys.exit(1)
    
    if args.position_threshold is not None:
        print(f"Position threshold: {args.position_threshold}")
    
    print("-" * 60)
    
    # Process structures
    all_results = process_structures(args.str_file, args.ss_folder, args.dots, args.dots_only, positions, args.position_threshold)
    
    if not all_results:
        print("No .ss_detected files found.")
        sys.exit(1)
    
    total_files = sum(len(results) for results in all_results.values())
    print("-" * 60)
    print(f"Processed {total_files} structures across {len(all_results)} trajectories.")
    
    # Save CSV
    if args.csv:
        fieldnames = ['trajectory', 'file', 'frame_number', 'total_score', 'total_possible', 
                     'pair_score', 'dot_score', 'total_pairs', 'total_dots', 
                     'accuracy', 'pair_accuracy', 'dot_accuracy']
        
        # Add position-specific fields if positions were specified
        if positions:
            fieldnames.extend(['position_accuracy', 'position_correct', 'position_total'])
        
        with open(args.csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for suffix, results in all_results.items():
                for result in results:
                    row = result.copy()
                    row['trajectory'] = suffix
                    writer.writerow(row)
        
        print(f"Results saved to {args.csv}")
    
    # Create plots
    if args.all:
        created_plot = plot_all_results_combined(all_results, args.output, args.dots, args.dots_only, positions_only)
        print(f"\nCreated combined plot: {created_plot}")

if __name__ == "__main__":
    main()