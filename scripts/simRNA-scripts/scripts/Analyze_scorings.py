#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import glob
import shutil

def parse_positions(position_spec):
    """
    Parse position specification from string.
    Supports formats like:
    - "24-29" (range)
    - "24,25,26,27" (comma-separated)
    - "24-29,35,40-42" (mixed)
    Returns list of 0-based indices
    """
    positions = []
    if not position_spec:
        return positions
    
    parts = position_spec.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            positions.extend(range(start-1, end))  # Convert to 0-based
        else:
            positions.append(int(part) - 1)  # Convert to 0-based
    
    return sorted(list(set(positions)))

def get_position_suffix(positions):
    """Generate a suffix string for position-specific files"""
    if not positions:
        return ""
    
    # Convert back to 1-based for the suffix
    pos_1based = [p + 1 for p in positions]
    
    # Create a compact representation
    ranges = []
    start = pos_1based[0]
    end = start
    
    for i in range(1, len(pos_1based)):
        if pos_1based[i] == end + 1:
            end = pos_1based[i]
        else:
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{end}")
            start = end = pos_1based[i]
    
    # Add the last range
    if start == end:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}-{end}")
    
    return "_pos_" + "_".join(ranges).replace("-", "to")

def run_command(cmd, quiet=False):
    """Run a command and return success status"""
    if not quiet:
        print(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {' '.join(cmd)}")
        if e.stdout:
            print(f"STDOUT: {e.stdout}")
        if e.stderr:
            print(f"STDERR: {e.stderr}")
        print(f"Return code: {e.returncode}")
        return False

def cleanup_pdbs():
    """Remove pdbs directory if it exists"""
    if os.path.exists("pdbs"):
        shutil.rmtree("pdbs")

def run_connectivity_analysis_for_temp(name, temp_level, pos_suffix="", quiet=False):
    """Run strand connectivity analysis for a single temperature level"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    connectivity_script = os.path.join(script_dir, "analyze_strand_connectivity.py")
    
    if not os.path.exists(connectivity_script):
        if not quiet:
            print(f"Warning: {connectivity_script} not found, skipping connectivity analysis")
        return None
    
    cmd = [
        "python", connectivity_script,
        "--name", name,
        "--min-temp", str(temp_level),
        "--max-temp", str(temp_level),
        "--positions-suffix", pos_suffix,
        "--skip-plots", # Don't create plots for individual temps
        "--debug"
    ]
    
    if quiet:
        cmd.append("--quiet")
    
    success = run_command(cmd, quiet)
    
    # Read and return the connectivity data
    if success:
        csv_file = f"{name}_connectivity_summary{pos_suffix}.csv"
        if os.path.exists(csv_file):
            try:
                df = pd.read_csv(csv_file)
                if len(df) > 0:
                    return df.iloc[0].to_dict()
            except Exception as e:
                if not quiet:
                    print(f"Error reading connectivity summary: {e}")
    
    return None

def process_temperature_level(temp_level, name, structure_pattern, min_frame, script_dir, 
                            positions=None, position_threshold=None, 
                            skip_connectivity=False, quiet=False):
    """Process a single temperature level with all three scoring methods and connectivity analysis"""
    if not quiet:
        print(f"Processing temperature level {temp_level}...")
    
    cleanup_pdbs()
    
    # Get position suffix for file naming
    pos_suffix = get_position_suffix(positions)
    
    # Extract temperature frames
    trafl_output = f"level{temp_level}_temp_frames{pos_suffix}.trafl"
    extract_script = os.path.join(script_dir, "extract_low_temp_frames.py")
    extract_cmd = [
        "python", extract_script, f"{name}.log",
        "--base-name", name, "--min-temp", str(temp_level),
        "--max-temp", str(temp_level), "--min-frame", str(min_frame),
        "--output", trafl_output
    ]
    
    if not run_command(extract_cmd, quiet) or not os.path.exists(trafl_output):
        return False, None
    
    # Extract PDB files (this creates .ss_detected files)
    extract_pdbs_script = os.path.join(script_dir, "extract_all_pdbs.py")
    extract_pdbs_cmd = [
        "python", extract_pdbs_script, "./",
        "--name", name, "--trafl", trafl_output
    ]
    
    if not run_command(extract_pdbs_cmd, quiet):
        return False, None
    
    # Compare structures for all three scoring methods
    compare_script = os.path.join(script_dir, "compare_ss.py")
    scoring_methods = [
        ("pairs", []),
        ("dots", ["--dots"]),
        ("dots_only", ["--dots-only"])
    ]
    
    for method, extra_args in scoring_methods:
        output_prefix = f"TEMP{temp_level}_{method}{pos_suffix}"
        csv_filename = f"accuracy_TEMP{temp_level}_{method}{pos_suffix}.csv"
        
        compare_cmd = [
            "python", compare_script, structure_pattern, "pdbs",
            "--output", output_prefix,
            "--csv", csv_filename,
            "--all"
        ] + extra_args
        
        # Add position-specific arguments
        if positions:
            compare_cmd.extend(["--positions", ",".join(map(lambda x: str(x+1), positions))])  # Convert back to 1-based
        if position_threshold is not None:
            compare_cmd.extend(["--position-threshold", str(position_threshold)])
        
        if not run_command(compare_cmd, quiet):
            return False, None
    
    # Run connectivity analysis BEFORE cleaning up pdbs
    connectivity_data = None
    if not skip_connectivity:
        if not quiet:
            print(f"  Running connectivity analysis for temperature {temp_level}...")
        connectivity_data = run_connectivity_analysis_for_temp(name, temp_level, pos_suffix, quiet)
    
    cleanup_pdbs()
    return True, connectivity_data

def calculate_stats(successful_levels, scoring_type, threshold=0.8, positions_only=False, pos_suffix=""):
    """Calculate accuracy statistics for temperature levels"""
    all_accuracies = []
    level_stats = {}
    
    for level in successful_levels:
        csv_file = f"accuracy_TEMP{level}_{scoring_type}{pos_suffix}.csv"
        if not os.path.exists(csv_file):
            continue
            
        try:
            df = pd.read_csv(csv_file)
            # Choose the appropriate accuracy column
            if positions_only and 'position_accuracy' in df.columns:
                accuracy_col = 'position_accuracy'
            else:
                accuracy_col = 'accuracy'
            
            if accuracy_col not in df.columns:
                continue
                
            accuracies = df[accuracy_col].values
            all_accuracies.extend(accuracies)
            
            total_frames = len(accuracies)
            frames_above_threshold = sum(acc >= threshold for acc in accuracies)
            percentage_above = (frames_above_threshold / total_frames * 100) if total_frames > 0 else 0
            
            level_stats[level] = {
                'percentage_above': percentage_above,
                'mean_accuracy': np.mean(accuracies),
                'median_accuracy': np.median(accuracies)
            }
        except Exception as e:
            print(f"Error processing CSV {csv_file}: {e}")
            continue
    
    # Overall stats
    if all_accuracies:
        total_all = len(all_accuracies)
        above_threshold = sum(acc >= threshold for acc in all_accuracies)
        overall_stats = {
            'percentage_above': (above_threshold / total_all * 100) if total_all > 0 else 0,
            'mean_accuracy': np.mean(all_accuracies),
            'median_accuracy': np.median(all_accuracies)
        }
    else:
        overall_stats = {'percentage_above': 0, 'mean_accuracy': 0, 'median_accuracy': 0}
    
    return overall_stats, level_stats

def update_summary_file(summary_file, label, successful_levels, threshold=0.8, use_median=False, 
                       use_mean=False, positions_only=False, pos_suffix=""):
    """Update summary CSV with results from all three scoring methods"""
    # Calculate stats for all three methods
    pairs_stats, pairs_level = calculate_stats(successful_levels, "pairs", threshold, positions_only, pos_suffix)
    dots_stats, dots_level = calculate_stats(successful_levels, "dots", threshold, positions_only, pos_suffix)
    dots_only_stats, dots_only_level = calculate_stats(successful_levels, "dots_only", threshold, positions_only, pos_suffix)
    
    # Create rows for each method
    rows = []
    for method, level_stats in [("pairs", pairs_level), ("dots", dots_level), ("dots_only", dots_only_level)]:
        row = {'label': f"{label}_{method}{pos_suffix}"}
        for level in successful_levels:
            if level in level_stats:
                if use_median:
                    row[str(level)] = round(level_stats[level]['median_accuracy'], 4)
                elif use_mean:
                    row[str(level)] = round(level_stats[level]['mean_accuracy'], 4)
                else:
                    row[str(level)] = round(level_stats[level]['percentage_above'], 1)
            else:
                row[str(level)] = None
        rows.append(row)
    
    # Read existing data or create new
    if os.path.exists(summary_file):
        try:
            existing_df = pd.read_csv(summary_file)
            new_df = pd.DataFrame(rows)
            updated_df = pd.concat([existing_df, new_df], ignore_index=True)
            
            # Get all temperature columns
            temp_cols = set()
            for col in updated_df.columns:
                if col != 'label' and col.isdigit():
                    temp_cols.add(col)
            
            final_columns = ['label'] + sorted(temp_cols, key=int)
            for col in final_columns:
                if col not in updated_df.columns:
                    updated_df[col] = None
            updated_df = updated_df[final_columns]
        except Exception as e:
            print(f"Error updating summary file: {e}")
            updated_df = pd.DataFrame(rows)
    else:
        updated_df = pd.DataFrame(rows)
    
    updated_df.to_csv(summary_file, index=False)
    print(f"Summary updated: {summary_file}")

def create_violin_plot(successful_levels, name, scoring_type="pairs", positions_only=False, pos_suffix=""):
    """Create violin plot for a scoring method"""
    all_data = []
    
    for level in successful_levels:
        csv_file = f"accuracy_TEMP{level}_{scoring_type}{pos_suffix}.csv"
        if not os.path.exists(csv_file):
            continue
            
        try:
            df = pd.read_csv(csv_file)
            # Choose the appropriate accuracy column
            if positions_only and 'position_accuracy' in df.columns:
                accuracy_col = 'position_accuracy'
            else:
                accuracy_col = 'accuracy'
            
            if accuracy_col in df.columns:
                df['Temperature'] = level
                df['plot_accuracy'] = df[accuracy_col]
                all_data.append(df)
        except Exception as e:
            print(f"Error processing {csv_file} for violin plot: {e}")
            continue
    
    if not all_data:
        return None
    
    combined_df = pd.concat(all_data, ignore_index=True)
    
    plt.figure(figsize=(12, 8))
    sns.violinplot(data=combined_df, x='Temperature', y='plot_accuracy', inner='box')
    
    scoring_titles = {"pairs": "Base Pairs Only", "dots": "Including Dots", "dots_only": "Dots Only"}
    title = scoring_titles.get(scoring_type, scoring_type.title())
    pos_title_suffix = " (Position-Specific)" if positions_only else ""
    
    plt.title(f'Accuracy Distribution - {title}{pos_title_suffix}\n({name})', fontsize=14, fontweight='bold')
    plt.xlabel('Temperature Level', fontsize=12)
    plt.ylabel('Accuracy', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.ylim(-0.05, 1.05)
    
    output_file = f"{name}_temperature_violin_{scoring_type}{pos_suffix}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_file

def create_comparison_plot(successful_levels, name, positions_only=False, pos_suffix=""):
    """Create comparison plot with all three scoring methods"""
    all_data = []
    
    for scoring_type in ["pairs", "dots", "dots_only"]:
        for level in successful_levels:
            csv_file = f"accuracy_TEMP{level}_{scoring_type}{pos_suffix}.csv"
            if not os.path.exists(csv_file):
                continue
                
            try:
                df = pd.read_csv(csv_file)
                # Choose the appropriate accuracy column
                if positions_only and 'position_accuracy' in df.columns:
                    accuracy_col = 'position_accuracy'
                else:
                    accuracy_col = 'accuracy'
                
                if accuracy_col in df.columns:
                    df['Temperature'] = level
                    scoring_map = {"pairs": "Base Pairs Only", "dots": "Including Dots", "dots_only": "Dots Only"}
                    df['Scoring'] = scoring_map[scoring_type]
                    df['plot_accuracy'] = df[accuracy_col]
                    all_data.append(df)
            except Exception as e:
                print(f"Error processing {csv_file} for comparison plot: {e}")
                continue
    
    if not all_data:
        return None
    
    combined_df = pd.concat(all_data, ignore_index=True)
    
    plt.figure(figsize=(16, 8))
    sns.violinplot(data=combined_df, x='Temperature', y='plot_accuracy', hue='Scoring', 
                   split=False, inner='box', palette=['steelblue', 'orange', 'green'])
    
    pos_title_suffix = " (Position-Specific)" if positions_only else ""
    plt.title(f'Accuracy Comparison: All Scoring Methods{pos_title_suffix}\n({name})', fontsize=14, fontweight='bold')
    plt.xlabel('Temperature Level', fontsize=12)
    plt.ylabel('Accuracy', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.ylim(-0.05, 1.05)
    plt.legend(title='Scoring Method', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    output_file = f"{name}_temperature_comparison{pos_suffix}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_file

def create_connectivity_plots(connectivity_data, name, pos_suffix=""):
    """Create connectivity plots from collected data"""
    if not connectivity_data:
        return []
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Import the plotting function from analyze_strand_connectivity
    sys.path.insert(0, script_dir)
    try:
        from analyze_strand_connectivity import create_connectivity_plots as make_plots
        from analyze_strand_connectivity import save_connectivity_summary
        
        # Save the combined summary
        save_connectivity_summary(connectivity_data, name, pos_suffix)
        
        # Create plots
        return make_plots(connectivity_data, name, pos_suffix)
    except ImportError as e:
        print(f"Warning: Could not import connectivity plotting functions: {e}")
        return []
    finally:
        sys.path.pop(0)

def create_pdf_report(name, positions_only=False, pos_suffix=""):
    """Create combined PDF report"""
    png_files = [f for f in glob.glob("*.png") if any(x in f for x in ["TEMP", name]) and pos_suffix in f]
    
    if not png_files:
        return None
    
    # Sort files by temperature and type
    def sort_key(filename):
        import re
        temp_match = re.search(r'TEMP(\d+)', filename)
        if temp_match:
            temp_num = int(temp_match.group(1))
            if 'pairs' in filename:
                return (0, temp_num, 0)
            elif 'dots' in filename:
                return (0, temp_num, 1)
            else:
                return (0, temp_num, 2)
        elif 'connectivity' in filename:
            return (1, 0, 0)
        elif 'violin' in filename:
            return (2, 0, 0)
        elif 'comparison' in filename:
            return (3, 0, 0)
        else:
            return (4, 0, 0)
    
    png_files.sort(key=sort_key)
    
    pdf_filename = f"{name}_analysis_report{pos_suffix}.pdf"
    
    with PdfPages(pdf_filename) as pdf:
        for png_file in png_files:
            if os.path.exists(png_file):
                img = plt.imread(png_file)
                fig, ax = plt.subplots(figsize=(11, 8.5))
                ax.imshow(img)
                ax.axis('off')
                fig.suptitle(os.path.basename(png_file), fontsize=12, fontweight='bold')
                pdf.savefig(fig, bbox_inches='tight', dpi=150)
                plt.close(fig)
    
    return pdf_filename

def main():
    parser = argparse.ArgumentParser(description='RNA structure temperature analysis with position-specific options')
    parser.add_argument('--name', '-n', required=True, help='Base name for files')
    parser.add_argument('--min-temp', '-t', type=int, required=True, help='Min temperature level')
    parser.add_argument('--max-temp', '-T', type=int, required=True, help='Max temperature level')
    parser.add_argument('--structure', '-s', required=True, help='Structure file pattern')
    parser.add_argument('--min-frame', '-f', type=int, default=1000, help='Min frame number (default: 1000)')
    parser.add_argument('--positions', '-p', help='Positions to analyze (e.g., "24-29" or "24,25,26,27")')
    parser.add_argument('--position-threshold', type=float, 
                       help='Minimum position accuracy threshold to include frames')
    parser.add_argument('--continue-on-error', action='store_true', help='Continue if one level fails')
    parser.add_argument('--skip-plots', action='store_true', help='Skip generating plots')
    parser.add_argument('--skip-connectivity', action='store_true', help='Skip strand connectivity analysis')
    parser.add_argument('--summary', help='Summary CSV file to update')
    parser.add_argument('--label', help='Label for summary file')
    parser.add_argument('--threshold', type=float, default=0.8, help='Accuracy threshold (default: 0.8)')
    parser.add_argument('--median', action='store_true', help='Use median in summary')
    parser.add_argument('--mean', action='store_true', help='Use mean in summary')
    parser.add_argument('--quiet', '-q', action='store_true', help='Reduce output')
    
    args = parser.parse_args()
    
    # Basic validation
    if args.summary and not args.label:
        print("Error: --label required with --summary")
        sys.exit(1)
    
    if args.min_temp > args.max_temp:
        print("Error: min-temp cannot be greater than max-temp")
        sys.exit(1)
    
    # Parse positions
    positions = parse_positions(args.positions) if args.positions else None
    positions_only = positions is not None
    pos_suffix = get_position_suffix(positions)
    
    if positions and not args.quiet:
        print(f"Analyzing positions: {[p+1 for p in positions]} (1-based)")
        print(f"Position suffix: {pos_suffix}")
    
    if args.position_threshold is not None and not args.quiet:
        print(f"Position threshold: {args.position_threshold}")
    
    # Check required files
    log_file = f"{args.name}.log"
    if not os.path.exists(log_file):
        print(f"Error: {log_file} not found")
        sys.exit(1)
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    required_scripts = ["extract_low_temp_frames.py", "extract_all_pdbs.py", "compare_ss.py"]
    for script in required_scripts:
        if not os.path.exists(os.path.join(script_dir, script)):
            print(f"Error: {script} not found in {script_dir}")
            sys.exit(1)
    
    if not args.quiet:
        print(f"Analyzing {args.name}: levels {args.min_temp}-{args.max_temp}")
    
    # Process temperature levels
    successful_levels = []
    failed_levels = []
    connectivity_data = []
    
    for temp_level in range(args.min_temp, args.max_temp + 1):
        success, conn_data = process_temperature_level(
            temp_level, args.name, args.structure, 
            args.min_frame, script_dir, positions, 
            args.position_threshold, args.skip_connectivity, args.quiet
        )
        
        if success:
            successful_levels.append(temp_level)
            if conn_data:
                connectivity_data.append(conn_data)
        else:
            failed_levels.append(temp_level)
            if not args.continue_on_error:
                break
    
    cleanup_pdbs()
    
    print(f"Completed: {len(successful_levels)}/{args.max_temp - args.min_temp + 1} levels")
    if failed_levels:
        print(f"Failed levels: {failed_levels}")
    
    # Update summary file
    if args.summary and successful_levels:
        update_summary_file(args.summary, args.label, successful_levels, 
                          args.threshold, args.median, args.mean, positions_only, pos_suffix)
    
    # Create connectivity plots from collected data
    if not args.skip_connectivity and connectivity_data:
        if not args.quiet:
            print("\nCreating connectivity plots...")
        create_connectivity_plots(connectivity_data, args.name, pos_suffix)
    
    # Create plots
    if not args.skip_plots and successful_levels:
        for method in ["pairs", "dots", "dots_only"]:
            create_violin_plot(successful_levels, args.name, method, positions_only, pos_suffix)
        create_comparison_plot(successful_levels, args.name, positions_only, pos_suffix)
        pdf_file = create_pdf_report(args.name, positions_only, pos_suffix)
        if pdf_file:
            print(f"PDF report: {pdf_file}")
    
    if failed_levels:
        sys.exit(1)
    else:
        print("✓ Analysis complete")

if __name__ == "__main__":
    main()