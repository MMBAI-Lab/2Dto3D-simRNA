#!/usr/bin/env python3
import argparse
import sys
import os
import pandas as pd
import numpy as np
import glob
import re
import time
import random

def acquire_file_lock(file_path, timeout=30):
    """Acquire an exclusive lock on a file with timeout"""
    lock_file = f"{file_path}.lock"
    start_time = time.time()
    
    while time.time() - start_time < timeout:
        try:
            lock_fd = os.open(lock_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(lock_fd, str(os.getpid()).encode())
            return lock_fd, lock_file
        except FileExistsError:
            time.sleep(0.1 + random.uniform(0, 0.1))
            continue
    
    raise TimeoutError(f"Could not acquire lock on {file_path} within {timeout} seconds")

def release_file_lock(lock_fd, lock_file):
    """Release the file lock"""
    try:
        os.close(lock_fd)
        os.unlink(lock_file)
    except:
        pass

def calculate_accuracy_stats(successful_levels, scoring_type, threshold=0.8):
    """Calculate accuracy statistics for all temperature levels"""
    all_accuracies = []
    level_stats = {}
    
    for level in successful_levels:
        csv_file = f"accuracy_TEMP{level}_{scoring_type}.csv"
        
        if not os.path.exists(csv_file):
            continue
        
        try:
            df = pd.read_csv(csv_file)
            if 'accuracy' not in df.columns:
                continue
            
            accuracies = df['accuracy'].values
            all_accuracies.extend(accuracies)
            
            total_frames = len(accuracies)
            frames_above_threshold = sum(acc >= threshold for acc in accuracies)
            percentage_above = (frames_above_threshold / total_frames * 100) if total_frames > 0 else 0
            
            level_stats[level] = {
                'total_frames': total_frames,
                'frames_above_threshold': frames_above_threshold,
                'percentage_above': percentage_above,
                'mean_accuracy': np.mean(accuracies),
                'median_accuracy': np.median(accuracies),
                'std_accuracy': np.std(accuracies)
            }
        except:
            continue
    
    # Overall stats
    if all_accuracies:
        total_all_frames = len(all_accuracies)
        frames_above_threshold_all = sum(acc >= threshold for acc in all_accuracies)
        percentage_above_all = (frames_above_threshold_all / total_all_frames * 100) if total_all_frames > 0 else 0
        
        overall_stats = {
            'total_frames': total_all_frames,
            'frames_above_threshold': frames_above_threshold_all,
            'percentage_above': percentage_above_all,
            'mean_accuracy': np.mean(all_accuracies),
            'median_accuracy': np.median(all_accuracies),
            'std_accuracy': np.std(all_accuracies)
        }
    else:
        overall_stats = {
            'total_frames': 0,
            'frames_above_threshold': 0,
            'percentage_above': 0,
            'mean_accuracy': 0,
            'median_accuracy': 0,
            'std_accuracy': 0
        }
    
    return overall_stats, level_stats

def update_summary_file(summary_file, label, successful_levels, threshold=0.8, use_median=False, use_mean=False):
    """Update the summary CSV file with results from ALL THREE scoring methods"""
    lock_fd = None
    lock_file = None
    
    try:
        lock_fd, lock_file = acquire_file_lock(summary_file, timeout=60)
        
        # Calculate stats for all three scoring methods
        pairs_stats, pairs_level_stats = calculate_accuracy_stats(successful_levels, "pairs", threshold)
        dots_stats, dots_level_stats = calculate_accuracy_stats(successful_levels, "dots", threshold)
        dots_only_stats, dots_only_level_stats = calculate_accuracy_stats(successful_levels, "dots_only", threshold)
        
        # Create base summary rows for all three scoring methods
        pairs_row = {'label': f"{label}_pairs"}
        dots_row = {'label': f"{label}_dots"}
        dots_only_row = {'label': f"{label}_dots_only"}
        
        # Add either percentage above threshold, median values, or mean values for each successful temperature level
        for level in successful_levels:
            if level in pairs_level_stats:
                if use_median:
                    pairs_row[str(level)] = round(pairs_level_stats[level]['median_accuracy'], 4)
                elif use_mean:
                    pairs_row[str(level)] = round(pairs_level_stats[level]['mean_accuracy'], 4)
                else:
                    pairs_row[str(level)] = round(pairs_level_stats[level]['percentage_above'], 1)
            else:
                pairs_row[str(level)] = None
            
            if level in dots_level_stats:
                if use_median:
                    dots_row[str(level)] = round(dots_level_stats[level]['median_accuracy'], 4)
                elif use_mean:
                    dots_row[str(level)] = round(dots_level_stats[level]['mean_accuracy'], 4)
                else:
                    dots_row[str(level)] = round(dots_level_stats[level]['percentage_above'], 1)
            else:
                dots_row[str(level)] = None
            
            if level in dots_only_level_stats:
                if use_median:
                    dots_only_row[str(level)] = round(dots_only_level_stats[level]['median_accuracy'], 4)
                elif use_mean:
                    dots_only_row[str(level)] = round(dots_only_level_stats[level]['mean_accuracy'], 4)
                else:
                    dots_only_row[str(level)] = round(dots_only_level_stats[level]['percentage_above'], 1)
            else:
                dots_only_row[str(level)] = None
        
        # Define column order (label + actual temperature levels used)
        columns = ['label'] + [str(level) for level in sorted(successful_levels)]
        
        # Read existing data if file exists
        if os.path.exists(summary_file):
            try:
                existing_df = pd.read_csv(summary_file)
                
                # Add new rows
                new_rows = [pairs_row, dots_row, dots_only_row]
                new_df = pd.DataFrame(new_rows)
                updated_df = pd.concat([existing_df, new_df], ignore_index=True)
                
                # Get all unique temperature level columns
                all_temp_columns = set()
                for col in existing_df.columns:
                    if col != 'label' and col.isdigit():
                        all_temp_columns.add(col)
                for col in new_df.columns:
                    if col != 'label' and col.isdigit():
                        all_temp_columns.add(col)
                
                # Create final column order
                final_columns = ['label'] + sorted(all_temp_columns, key=int)
                
                # Ensure all columns are present
                for col in final_columns:
                    if col not in updated_df.columns:
                        updated_df[col] = None
                
                updated_df = updated_df[final_columns]
                
            except:
                new_rows = [pairs_row, dots_row, dots_only_row]
                updated_df = pd.DataFrame(new_rows)
                updated_df = updated_df.reindex(columns=columns)
        else:
            # Create new DataFrame
            new_rows = [pairs_row, dots_row, dots_only_row]
            updated_df = pd.DataFrame(new_rows)
            updated_df = updated_df.reindex(columns=columns)
        
        # Save updated summary
        updated_df.to_csv(summary_file, index=False)
        
        # Print summary stats
        if use_median:
            metric_name = "median accuracy"
            print(f"Overall median accuracy:")
            print(f"  Pairs only:  {pairs_stats['median_accuracy']:.4f}")
            print(f"  With dots:   {dots_stats['median_accuracy']:.4f}")
            print(f"  Dots only:   {dots_only_stats['median_accuracy']:.4f}")
        elif use_mean:
            metric_name = "mean accuracy"
            print(f"Overall mean accuracy:")
            print(f"  Pairs only:  {pairs_stats['mean_accuracy']:.4f}")
            print(f"  With dots:   {dots_stats['mean_accuracy']:.4f}")
            print(f"  Dots only:   {dots_only_stats['mean_accuracy']:.4f}")
        else:
            metric_name = f"percentage above {threshold*100:.0f}% threshold"
            print(f"Overall stats:")
            print(f"  Pairs only:  {pairs_stats['percentage_above']:.1f}% ({pairs_stats['frames_above_threshold']}/{pairs_stats['total_frames']}) frames above threshold")
            print(f"  With dots:   {dots_stats['percentage_above']:.1f}% ({dots_stats['frames_above_threshold']}/{dots_stats['total_frames']}) frames above threshold")
            print(f"  Dots only:   {dots_only_stats['percentage_above']:.1f}% ({dots_only_stats['frames_above_threshold']}/{dots_only_stats['total_frames']}) frames above threshold")
        
        return True
        
    except Exception as e:
        print(f"Error updating summary file: {e}")
        return False
        
    finally:
        if lock_fd is not None and lock_file is not None:
            release_file_lock(lock_fd, lock_file)

def auto_detect_levels(directory="."):
    """Auto-detect temperature levels from existing CSV files"""
    pairs_files = glob.glob(os.path.join(directory, "accuracy_TEMP*_pairs.csv"))
    dots_files = glob.glob(os.path.join(directory, "accuracy_TEMP*_dots.csv"))
    dots_only_files = glob.glob(os.path.join(directory, "accuracy_TEMP*_dots_only.csv"))
    
    temp_levels = set()
    for filename in pairs_files + dots_files + dots_only_files:
        match = re.search(r'accuracy_TEMP(\d+)_', os.path.basename(filename))
        if match:
            temp_levels.add(int(match.group(1)))
    
    return sorted(temp_levels)

def main():
    parser = argparse.ArgumentParser(description='Create summary file from existing accuracy CSV files (all three scoring methods)')
    parser.add_argument('--summary', required=True, help='Summary CSV file to update with results')
    parser.add_argument('--label', required=True, help='Label for this run in the summary file')
    parser.add_argument('--levels', type=int, nargs='+', help='Temperature levels to process (e.g., 1 2 3 4 5). If not specified, will auto-detect.')
    parser.add_argument('--threshold', type=float, default=0.8, help='Accuracy threshold for summary statistics (default: 0.8 = 80%%)')
    parser.add_argument('--median', action='store_true', help='Store median accuracy values in summary instead of percentage above threshold')
    parser.add_argument('--mean', action='store_true', help='Store mean accuracy values in summary instead of percentage above threshold')
    parser.add_argument('--directory', default='.', help='Directory containing accuracy CSV files (default: current directory)')
    
    args = parser.parse_args()
    
    # Basic validation
    if args.threshold < 0 or args.threshold > 1:
        print("Error: threshold must be between 0 and 1")
        sys.exit(1)
    
    if args.median and args.mean:
        print("Error: --median and --mean cannot be used together")
        sys.exit(1)
    
    # Change to specified directory
    if args.directory != '.':
        if not os.path.exists(args.directory):
            print(f"Error: Directory {args.directory} does not exist")
            sys.exit(1)
        original_dir = os.getcwd()
        os.chdir(args.directory)
    
    # Determine temperature levels
    if args.levels:
        successful_levels = sorted(args.levels)
    else:
        successful_levels = auto_detect_levels()
        if not successful_levels:
            print("Error: No temperature levels found. Please specify --levels manually.")
            sys.exit(1)
    
    # Check if required CSV files exist
    missing_files = []
    for level in successful_levels:
        for scoring_type in ["pairs", "dots", "dots_only"]:
            csv_file = f"accuracy_TEMP{level}_{scoring_type}.csv"
            if not os.path.exists(csv_file):
                missing_files.append(csv_file)
    
    if missing_files:
        print("Error: Missing required CSV files:")
        for file in missing_files:
            print(f"  {file}")
        sys.exit(1)
    
    print(f"Creating Summary for Label: {args.label}")
    print(f"Temperature levels: {successful_levels}")
    if args.median:
        print(f"Summary metric: median accuracy")
    elif args.mean:
        print(f"Summary metric: mean accuracy")
    else:
        print(f"Summary metric: percentage above {args.threshold} ({args.threshold*100:.1f}%) threshold")
    print(f"Summary file: {args.summary}")
    
    # Update summary file
    success = update_summary_file(args.summary, args.label, successful_levels, args.threshold, args.median, args.mean)
    
    if success:
        print(f"\n✓ Summary file updated successfully!")
    else:
        print(f"\n✗ Failed to update summary file")
        sys.exit(1)
    
    # Change back to original directory if needed
    if args.directory != '.':
        os.chdir(original_dir)

if __name__ == "__main__":
    main()