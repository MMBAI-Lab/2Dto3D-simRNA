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
import re
from collections import defaultdict

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

def normalize_structure(structure):
    """
    Normalize a structure to a canonical single-line format.
    Handles multi-line structures by joining lines and removing whitespace.
    Now handles extended bracket notation: ()[]{}
    """
    if pd.isna(structure) or structure == '':
        return ''
    
    # Convert to string and normalize
    structure_str = str(structure)
    
    # Remove whitespace and newlines, keep structure characters including extended brackets and '+' for spaces
    normalized = ''.join(char for char in structure_str if char in '()[]{}.-+<>')
    
    return normalized

def calculate_structure_similarity(struct1, struct2):
    """
    Calculate similarity between two structures.
    Returns similarity as a fraction (0.0 to 1.0).
    """
    if not struct1 or not struct2:
        return 0.0
    
    if len(struct1) != len(struct2):
        return 0.0
    
    matches = sum(1 for a, b in zip(struct1, struct2) if a == b)
    return matches / len(struct1)

def collect_structures_from_csvs(successful_levels, scoring_methods, pos_suffix="", target_structure=None):
    """
    Collect all structures from generated CSV files and create unique structure summaries.
    """
    all_structures = []
    structure_metadata = defaultdict(lambda: {
        'count': 0,
        'temperature_levels': set(),
        'scoring_methods': set(),
        'accuracies': [],
        'first_seen_temp': None,
        'first_seen_method': None
    })
    
    # Normalize target structure if provided
    normalized_target = normalize_structure(target_structure) if target_structure else None
    
    for level in successful_levels:
        for method in scoring_methods:
            csv_file = f"accuracy_TEMP{level}_{method}{pos_suffix}.csv"
            if not os.path.exists(csv_file):
                print(f"Warning: CSV file {csv_file} not found")
                continue
                
            try:
                df = pd.read_csv(csv_file)
                print(f"Processing {csv_file} with columns: {list(df.columns)}")  # Debug
                
                # Look for structure columns with more specific patterns
                structure_cols = []
                
                # First, look for columns with obvious structure names
                for col in df.columns:
                    col_lower = col.lower()
                    if any(keyword in col_lower for keyword in ['structure', 'predicted', 'simulated', 'secondary', 'ss_']):
                        structure_cols.append(col)
                        print(f"Found structure column by name: {col}")  # Debug
                
                # If no obvious structure column found, examine column content
                if not structure_cols:
                    print("No obvious structure columns found, examining content...")  # Debug
                    for col in df.columns:
                        if col.lower() in ['accuracy', 'position_accuracy', 'frame', 'temperature', 'level']:
                            continue
                            
                        # Check if column contains structure-like data
                        sample_values = df[col].dropna().head(5)
                        if len(sample_values) > 0:
                            for sample_val in sample_values:
                                sample_str = str(sample_val)
                                # Check if it looks like a structure (has parentheses and reasonable length)
                                if ('(' in sample_str and ')' in sample_str and 
                                    len(sample_str) > 10 and 
                                    sum(1 for c in sample_str if c in '().-+') / len(sample_str) > 0.5):
                                    structure_cols.append(col)
                                    print(f"Found structure column by content: {col}, sample: {sample_str[:50]}")  # Debug
                                    break
                
                if not structure_cols:
                    print(f"No structure columns found in {csv_file}")
                    print(f"Available columns: {list(df.columns)}")
                    print(f"Sample data from first few columns:")
                    for col in df.columns[:5]:
                        sample = df[col].dropna().iloc[0] if len(df[col].dropna()) > 0 else 'N/A'
                        print(f"  {col}: {sample}")
                    continue
                
                # Process each structure column found
                for struct_col in structure_cols:
                    print(f"Processing structure column: {struct_col}")  # Debug
                    valid_structures = 0
                    
                    for idx, row in df.iterrows():
                        structure = row.get(struct_col, '')
                        accuracy = row.get('accuracy', row.get('position_accuracy', 0))
                        
                        normalized_struct = normalize_structure(structure)
                        
                        # Skip empty or too short structures
                        if not normalized_struct or len(normalized_struct) < 5:
                            continue
                        
                        valid_structures += 1
                        if valid_structures <= 3:  # Debug first few structures
                            print(f"  Structure {valid_structures}: '{structure}' -> '{normalized_struct}'")
                        
                        # Calculate similarity to target if available
                        similarity_to_target = 0.0
                        matches_target = False
                        if normalized_target:
                            similarity_to_target = calculate_structure_similarity(normalized_struct, normalized_target)
                            matches_target = similarity_to_target >= 0.95  # 95% similarity threshold
                        
                        # Store individual structure record
                        structure_record = {
                            'structure': structure,
                            'normalized_structure': normalized_struct,
                            'temperature': level,
                            'scoring_method': method,
                            'accuracy': accuracy,
                            'similarity_to_target': similarity_to_target,
                            'matches_target': matches_target,
                            'source_file': csv_file,
                            'structure_column': struct_col
                        }
                        all_structures.append(structure_record)
                        
                        # Update metadata for unique structures
                        meta = structure_metadata[normalized_struct]
                        meta['count'] += 1
                        meta['temperature_levels'].add(level)
                        meta['scoring_methods'].add(method)
                        meta['accuracies'].append(accuracy)
                        
                        if meta['first_seen_temp'] is None or level < meta['first_seen_temp']:
                            meta['first_seen_temp'] = level
                            meta['first_seen_method'] = method

                    print(f"Found {valid_structures} valid structures in column {struct_col})")  # Debug

            except Exception as e:
                print(f"Error processing {csv_file} for structure collection: {e}")
                import traceback
                traceback.print_exc()
                continue
    
    print(f"Total structures collected: {len(all_structures)}")  # Debug
    print(f"Unique structures found: {len(structure_metadata)}")  # Debug
    
    return all_structures, structure_metadata

def create_unique_structures_summary(all_structures, structure_metadata, name, pos_suffix="", target_structure=None):
    """
    Create unique structures summary similar to nupack analyzer output.
    """
    if not all_structures:
        print("No structures found for summary")
        return None
    
    total_structures = len(all_structures)
    unique_structures_file = f"{name}_unique_structures{pos_suffix}.csv"
    all_structures_file = f"{name}_all_structures{pos_suffix}.csv"
    
    # Create unique structures summary
    unique_data = []
    for norm_struct, meta in structure_metadata.items():
        frequency = meta['count'] / total_structures
        avg_accuracy = np.mean(meta['accuracies']) if meta['accuracies'] else 0
        temp_range = f"{min(meta['temperature_levels'])}-{max(meta['temperature_levels'])}" if len(meta['temperature_levels']) > 1 else str(min(meta['temperature_levels']))
        methods = ','.join(sorted(meta['scoring_methods']))
        
        # Calculate similarity to target if available
        similarity_to_target = 0.0
        matches_target = False
        if target_structure:
            normalized_target = normalize_structure(target_structure)
            if normalized_target:
                similarity_to_target = calculate_structure_similarity(norm_struct, normalized_target)
                matches_target = similarity_to_target >= 0.95
        
        unique_data.append({
            'structure': norm_struct,
            'count': meta['count'],
            'frequency': frequency,
            'avg_accuracy': avg_accuracy,
            'temperature_range': temp_range,
            'scoring_methods': methods,
            'first_seen_temp': meta['first_seen_temp'],
            'first_seen_method': meta['first_seen_method'],
            'similarity_to_target': similarity_to_target,
            'matches_target': matches_target
        })
    
    # Sort by count (most frequent first)
    unique_data.sort(key=lambda x: x['count'], reverse=True)
    
    # Write unique structures CSV
    unique_df = pd.DataFrame(unique_data)
    unique_df.to_csv(unique_structures_file, index=False)
    
    # Write all structures CSV
    all_df = pd.DataFrame(all_structures)
    all_df.to_csv(all_structures_file, index=False)
    
    print(f"✓ Created {unique_structures_file} with {len(unique_data)} unique structures")
    print(f"✓ Created {all_structures_file} with {total_structures} total structure records")
    
    # Print summary statistics
    print(f"\nStructure Analysis Summary:")
    print(f"  Total structures analyzed: {total_structures}")
    print(f"  Unique structures found: {len(unique_data)}")
    print(f"  Structure diversity: {len(unique_data)/total_structures:.2%}")
    
    if target_structure:
        matching_count = sum(1 for s in unique_data if s['matches_target'])
        matching_frequency = sum(s['frequency'] for s in unique_data if s['matches_target'])
        print(f"  Structures matching target (≥95% similarity): {matching_count}")
        print(f"  Frequency of target-like structures: {matching_frequency:.2%}")
    
    # Top 5 most frequent structures
    print(f"\nTop 5 Most Frequent Structures:")
    for i, struct in enumerate(unique_data[:5]):
        print(f"  {i+1}. {struct['structure'][:50]}{'...' if len(struct['structure']) > 50 else ''}")
        print(f"     Count: {struct['count']} ({struct['frequency']:.2%}), Avg Accuracy: {struct['avg_accuracy']:.3f}")
        if target_structure:
            print(f"     Similarity to target: {struct['similarity_to_target']:.3f}")
    
    return unique_structures_file, all_structures_file

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

def process_temperature_level(temp_level, name, structure_pattern, min_frame, script_dir, 
                            positions=None, position_threshold=None, quiet=False,
                            collect_structures=False, structure_collector=None):
    """Process a single temperature level with all three scoring methods"""
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
        return False
    
    # Extract PDB files (this creates .ss_detected files)
    extract_pdbs_script = os.path.join(script_dir, "extract_all_pdbs.py")
    extract_pdbs_cmd = [
        "python", extract_pdbs_script, "./",
        "--name", name, "--trafl", trafl_output
    ]
    
    if not run_command(extract_pdbs_cmd, quiet):
        return False
    
    # Collect structures from .ss_detected files before they get deleted
    if collect_structures and structure_collector is not None:
        structure_collector.collect_from_pdbs_folder(temp_level, "pdbs")
    
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
            compare_cmd.extend(["--positions", ",".join(map(lambda x: str(x+1), positions))])
        if position_threshold is not None:
            compare_cmd.extend(["--position-threshold", str(position_threshold)])
        
        if not run_command(compare_cmd, quiet):
            return False
        
        # Add accuracy information to structure collector
        if collect_structures and structure_collector is not None:
            structure_collector.add_accuracy_data(temp_level, method, csv_filename)
    
    cleanup_pdbs()
    return True

class StructureCollector:
    """Collects structures during processing before files are cleaned up"""
    
    def __init__(self, target_structure=None):
        self.all_structures = []
        self.structure_metadata = defaultdict(lambda: {
            'count': 0,
            'temperature_levels': set(),
            'scoring_methods': set(),
            'accuracies': [],
            'first_seen_temp': None,
            'first_seen_method': None,
            'files': []
        })
        self.target_structure = target_structure
        self.normalized_target = normalize_structure(target_structure) if target_structure else None
        self.temp_structures = {}  # temp storage for structures by temperature
        
    def collect_from_pdbs_folder(self, temp_level, pdbs_folder):
        """Collect all structures from .ss_detected files in pdbs folder"""
        if not os.path.exists(pdbs_folder):
            return
            
        print(f"Collecting structures from {pdbs_folder} for temperature {temp_level}")
        temp_structures = {}
        
        for filename in os.listdir(pdbs_folder):
            if filename.endswith('.ss_detected'):
                file_path = os.path.join(pdbs_folder, filename)
                structure = read_structure_from_ss_detected(file_path)
                
                if structure:
                    normalized_struct = normalize_structure(structure)
                    if normalized_struct and len(normalized_struct) >= 5:
                        temp_structures[filename] = {
                            'structure': structure,
                            'normalized_structure': normalized_struct,
                            'temperature': temp_level,
                            'filename': filename
                        }
        
        self.temp_structures[temp_level] = temp_structures
        print(f"Collected {len(temp_structures)} structures for temperature {temp_level}")
    
    def add_accuracy_data(self, temp_level, method, csv_filename):
        """Add accuracy data from CSV to the collected structures"""
        if temp_level not in self.temp_structures:
            return
            
        if not os.path.exists(csv_filename):
            return
            
        try:
            df = pd.read_csv(csv_filename)
            if 'file' not in df.columns:
                return
                
            for idx, row in df.iterrows():
                ss_filename = row['file']
                accuracy = row.get('accuracy', row.get('position_accuracy', 0))
                
                if ss_filename in self.temp_structures[temp_level]:
                    struct_data = self.temp_structures[temp_level][ss_filename].copy()
                    struct_data['scoring_method'] = method
                    struct_data['accuracy'] = accuracy
                    struct_data['frame_number'] = row.get('frame_number', 0)
                    
                    # Calculate similarity to target
                    similarity_to_target = 0.0
                    matches_target = False
                    if self.normalized_target:
                        similarity_to_target = calculate_structure_similarity(
                            struct_data['normalized_structure'], self.normalized_target
                        )
                        matches_target = similarity_to_target >= 0.95
                    
                    struct_data['similarity_to_target'] = similarity_to_target
                    struct_data['matches_target'] = matches_target
                    
                    self.all_structures.append(struct_data)
                    
                    # Update metadata
                    normalized_struct = struct_data['normalized_structure']
                    meta = self.structure_metadata[normalized_struct]
                    meta['count'] += 1
                    meta['temperature_levels'].add(temp_level)
                    meta['scoring_methods'].add(method)
                    meta['accuracies'].append(accuracy)
                    meta['files'].append(ss_filename)
                    
                    if meta['first_seen_temp'] is None or temp_level < meta['first_seen_temp']:
                        meta['first_seen_temp'] = temp_level
                        meta['first_seen_method'] = method
                        
        except Exception as e:
            print(f"Error adding accuracy data from {csv_filename}: {e}")
    
    def get_results(self):
        """Return collected structures and metadata"""
        return self.all_structures, dict(self.structure_metadata)
    
def read_structure_from_ss_detected(file_path):
    """Read multi-line structure from .ss_detected file and return as extended dot-bracket notation"""
    try:
        with open(file_path, 'r') as f:
            structure_lines = []
            for line in f:
                line = line.strip()
                # Skip comment lines and empty lines
                if line and not line.startswith('//'):
                    # Replace spaces with '+' signs instead of removing them
                    clean_line = line.replace(' ', '+')
                    structure_lines.append(clean_line)
            
            if not structure_lines:
                return None
            
            # Convert multi-line to single-line with extended brackets
            return convert_multiline_to_extended_brackets(structure_lines)
            
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def convert_multiline_to_extended_brackets(structure_lines):
    """
    Convert multi-line structure to single-line with extended bracket notation.
    Line 1: () 
    Line 2: []
    Line 3: {}
    Line 4: <>
    etc.
    
    If multiple symbols compete for the same position, use the first one in line order.
    """
    if not structure_lines:
        return ""
    
    # Find the maximum length to ensure all lines are the same length
    max_length = max(len(line) for line in structure_lines)
    
    # Pad all lines to the same length with dots
    padded_lines = []
    for line in structure_lines:
        padded_line = line + '.' * (max_length - len(line))
        padded_lines.append(padded_line)
    
    # Define bracket pairs for each line
    bracket_pairs = [
        ('(', ')'),  # Line 1
        ('[', ']'),  # Line 2  
        ('{', '}'),  # Line 3
        ('<', '>'),  # Line 4
    ]
    
    # Initialize result structure with dots
    result = ['.'] * max_length
    
    # Process each line with its corresponding bracket style
    for line_idx, line in enumerate(padded_lines):
        if line_idx >= len(bracket_pairs):
            # If we have more than 4 lines, cycle through bracket types
            bracket_idx = line_idx % len(bracket_pairs)
        else:
            bracket_idx = line_idx
            
        open_bracket, close_bracket = bracket_pairs[bracket_idx]
        
        # Convert parentheses in this line to the appropriate bracket style
        for pos, char in enumerate(line):
            # Only overwrite if current position is a dot (first line takes precedence)
            if result[pos] == '.':
                if char == '(':
                    result[pos] = open_bracket
                elif char == ')':
                    result[pos] = close_bracket
                elif char in '.-+':  # Now includes '+' for spaces
                    result[pos] = char
                # If char is already some other bracket from a previous line, keep the original
    
    return ''.join(result)

def collect_structures_from_ss_detected_files(successful_levels, scoring_methods, pos_suffix="", target_structure=None):
    """
    Collect all structures from .ss_detected files created during processing.
    This mirrors how compare_ss.py reads the structures.
    """
    all_structures = []
    structure_metadata = defaultdict(lambda: {
        'count': 0,
        'temperature_levels': set(),
        'scoring_methods': set(),
        'accuracies': [],
        'first_seen_temp': None,
        'first_seen_method': None,
        'files': []
    })
    
    # Normalize target structure if provided
    normalized_target = normalize_structure(target_structure) if target_structure else None
    
    # For each temperature level, we need to look at the .ss_detected files that were processed
    for level in successful_levels:
        print(f"Collecting structures from temperature level {level}")
        
        # Look for accuracy CSV files to get the file lists that were processed
        for method in scoring_methods:
            csv_file = f"accuracy_TEMP{level}_{method}{pos_suffix}.csv"
            if not os.path.exists(csv_file):
                continue
                
            try:
                # Read the CSV to get the list of files that were processed
                df = pd.read_csv(csv_file)
                
                if 'file' not in df.columns:
                    print(f"Warning: No 'file' column found in {csv_file}")
                    continue
                
                print(f"Processing {len(df)} files from {csv_file}")
                
                for idx, row in df.iterrows():
                    ss_filename = row['file']
                    accuracy = row.get('accuracy', row.get('position_accuracy', 0))
                    
                    # Construct path to .ss_detected file (they were in pdbs folder during processing)
                    # But since we cleanup pdbs after each level, we need to figure out the original structure
                    # The CSV should have been created when the pdbs folder existed
                    
                    # For now, let's try to read from a pdbs folder if it exists, or skip
                    ss_file_path = os.path.join('pdbs', ss_filename)
                    if not os.path.exists(ss_file_path):
                        # Try without pdbs folder (in case files are in current directory)
                        ss_file_path = ss_filename
                        if not os.path.exists(ss_file_path):
                            continue
                    
                    # Read the structure from the .ss_detected file
                    structure = read_structure_from_ss_detected(ss_file_path)
                    if not structure:
                        continue
                    
                    normalized_struct = normalize_structure(structure)
                    if not normalized_struct or len(normalized_struct) < 5:
                        continue
                    
                    # Calculate similarity to target if available
                    similarity_to_target = 0.0
                    matches_target = False
                    if normalized_target:
                        similarity_to_target = calculate_structure_similarity(normalized_struct, normalized_target)
                        matches_target = similarity_to_target >= 0.95
                    
                    # Store individual structure record
                    structure_record = {
                        'structure': structure,
                        'normalized_structure': normalized_struct,
                        'temperature': level,
                        'scoring_method': method,
                        'accuracy': accuracy,
                        'similarity_to_target': similarity_to_target,
                        'matches_target': matches_target,
                        'source_file': ss_filename,
                        'frame_number': row.get('frame_number', 0)
                    }
                    all_structures.append(structure_record)
                    
                    # Update metadata for unique structures
                    meta = structure_metadata[normalized_struct]
                    meta['count'] += 1
                    meta['temperature_levels'].add(level)
                    meta['scoring_methods'].add(method)
                    meta['accuracies'].append(accuracy)
                    meta['files'].append(ss_filename)
                    
                    if meta['first_seen_temp'] is None or level < meta['first_seen_temp']:
                        meta['first_seen_temp'] = level
                        meta['first_seen_method'] = method
                        
            except Exception as e:
                print(f"Error processing {csv_file} for structure collection: {e}")
                import traceback
                traceback.print_exc()
                continue
    
    print(f"Total structures collected: {len(all_structures)}")
    print(f"Unique structures found: {len(structure_metadata)}")
    
    return all_structures, structure_metadata

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

def create_pdf_report(name, positions_only=False, pos_suffix=""):
    """Create combined PDF report"""
    png_files = [f for f in glob.glob("*.png") if any(x in f for x in ["TEMP", name]) and pos_suffix in f]
    
    if not png_files:
        return None
    
    # Sort files by temperature and type
    def sort_key(filename):
        temp_match = re.search(r'TEMP(\d+)', filename)
        if temp_match:
            temp_num = int(temp_match.group(1))
            if 'pairs' in filename:
                return (0, temp_num, 0)
            elif 'dots' in filename:
                return (0, temp_num, 1)
            else:
                return (0, temp_num, 2)
        elif 'violin' in filename:
            return (1, 0, 0)
        elif 'comparison' in filename:
            return (2, 0, 0)
        else:
            return (3, 0, 0)
    
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

from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform

class RNAStructureAnalyzer:
    def __init__(self, similarity_threshold=0.8, cluster_threshold=0.7):
        self.similarity_threshold = similarity_threshold
        self.cluster_threshold = cluster_threshold
        
    def base_pair_comparison_similarity(self, struct1, struct2):
        """Calculate similarity based on base pair comparison - treats all bracket types as equivalent"""
        def extract_base_pairs(structure):
            """Extract base pairs treating all bracket types as equivalent"""
            stack = []
            pairs = set()
            # All bracket types are treated as equivalent
            open_brackets = '([{<'
            close_brackets = ')]}>'
            bracket_map = {'(': ')', '[': ']', '{': '}', '<': '>'}
            
            for i, char in enumerate(structure):
                if char in open_brackets:
                    stack.append(i)
                elif char in close_brackets:
                    if stack:
                        j = stack.pop()
                        pairs.add((j, i))  # Just store positions, ignore bracket type
            return pairs
        
        def normalize_for_comparison(structure):
            """Normalize structure by converting all brackets to () for comparison"""
            normalized = ""
            for char in structure:
                if char in '([{<':
                    normalized += '('
                elif char in ')]}> ':
                    normalized += ')'
                else:
                    normalized += char
            return normalized
        
        # Handle length differences by padding shorter structure
        max_len = max(len(struct1), len(struct2))
        struct1_padded = struct1.ljust(max_len, '.')
        struct2_padded = struct2.ljust(max_len, '.')
        
        # Extract base pairs (treating all bracket types as equivalent)
        pairs1 = extract_base_pairs(struct1_padded)
        pairs2 = extract_base_pairs(struct2_padded)
        
        if len(pairs1) == 0 and len(pairs2) == 0:
            return 1.0
        if len(pairs1) == 0 or len(pairs2) == 0:
            return 0.0
        
        # Calculate Jaccard similarity based on pairing positions only
        intersection = len(pairs1.intersection(pairs2))
        union = len(pairs1.union(pairs2))
        
        jaccard_similarity = intersection / union if union > 0 else 0.0
        
        # Also calculate position-wise similarity after normalization
        norm1 = normalize_for_comparison(struct1_padded)
        norm2 = normalize_for_comparison(struct2_padded)
        position_matches = sum(1 for c1, c2 in zip(norm1, norm2) if c1 == c2)
        position_similarity = position_matches / max_len
        
        # Combine both similarities (base pair similarity is more important)
        combined_similarity = 0.8 * jaccard_similarity + 0.2 * position_similarity
        
        return combined_similarity

    def position_wise_similarity(self, struct1, struct2):
        """Calculate position-wise similarity treating all bracket types as equivalent"""
        def normalize_for_comparison(structure):
            """Normalize structure by converting all brackets to () for comparison"""
            normalized = ""
            for char in structure:
                if char in '([{<':
                    normalized += '('
                elif char in ')]}> ':
                    normalized += ')'
                else:
                    normalized += char
            return normalized
        
        max_len = max(len(struct1), len(struct2))
        struct1_padded = struct1.ljust(max_len, '.')
        struct2_padded = struct2.ljust(max_len, '.')
        
        # Normalize both structures for comparison
        norm1 = normalize_for_comparison(struct1_padded)
        norm2 = normalize_for_comparison(struct2_padded)
        
        matches = sum(1 for c1, c2 in zip(norm1, norm2) if c1 == c2)
        return matches / max_len

    def find_target_similarities(self, df, target_structure):
        """Find similarities to a target structure treating all bracket types as equivalent"""
        structures = df['structure'].tolist()
        similarities = []
        
        for struct in structures:
            bp_sim = self.base_pair_comparison_similarity(target_structure, struct)
            pos_sim = self.position_wise_similarity(target_structure, struct)
            combined_sim = 0.7 * bp_sim + 0.3 * pos_sim
            similarities.append({
                'base_pair_similarity': bp_sim,
                'position_similarity': pos_sim,
                'combined_similarity': combined_sim
            })
        
        # Add similarity scores to dataframe
        result_df = df.copy()
        result_df['bp_similarity_to_target'] = [s['base_pair_similarity'] for s in similarities]
        result_df['pos_similarity_to_target'] = [s['position_similarity'] for s in similarities]
        result_df['combined_similarity_to_target'] = [s['combined_similarity'] for s in similarities]
        
        # Sort by combined similarity
        result_df = result_df.sort_values('combined_similarity_to_target', ascending=False)
        
        return result_df

    def create_cluster_summary(self, df_clustered):
        """Create summary statistics for each cluster"""
        cluster_summary = []
        
        for cluster_id in sorted(df_clustered['cluster'].unique()):
            cluster_data = df_clustered[df_clustered['cluster'] == cluster_id]
            
            # Find representative structure (most common)
            representative_idx = cluster_data['count'].idxmax()
            representative_structure = cluster_data.loc[representative_idx, 'structure']
            
            # Calculate cluster statistics
            summary = {
                'cluster_id': cluster_id,
                'size': len(cluster_data),
                'total_count': cluster_data['count'].sum(),
                'avg_count': cluster_data['count'].mean(),
                'representative_structure': representative_structure,
                'avg_accuracy': cluster_data['avg_accuracy'].mean(),
                'structures_in_cluster': cluster_data['structure'].tolist(),
                'structure_diversity': self._calculate_cluster_diversity(cluster_data['structure'].tolist())
            }
            cluster_summary.append(summary)
        
        # Sort by total_count in descending order
        cluster_summary_df = pd.DataFrame(cluster_summary)
        cluster_summary_df = cluster_summary_df.sort_values('total_count', ascending=False)
        
        return cluster_summary_df

    def _calculate_cluster_diversity(self, structures):
        """Calculate diversity within a cluster"""
        if len(structures) <= 1:
            return 0.0
        
        total_similarity = 0
        comparisons = 0
        
        for i in range(len(structures)):
            for j in range(i + 1, len(structures)):
                similarity = self.base_pair_comparison_similarity(structures[i], structures[j])
                total_similarity += similarity
                comparisons += 1
        
        return total_similarity / comparisons if comparisons > 0 else 0.0

    def perform_clustering_analysis(self, unique_structures_file, name, pos_suffix=""):
        """Perform clustering analysis on unique structures"""
        try:
            df = pd.read_csv(unique_structures_file)
            structures = df['structure'].tolist()
            
            if len(structures) < 2:
                print("Not enough structures for clustering")
                return None
                
            print(f"Performing clustering analysis on {len(structures)} unique structures...")
            print("Note: Clustering treats all bracket types ([{< as equivalent to ( for similarity")
            
            # Calculate similarity matrix (using bracket-normalized comparison)
            n = len(structures)
            similarity_matrix = np.zeros((n, n))
            
            for i in range(n):
                for j in range(i, n):
                    if i == j:
                        similarity_matrix[i][j] = 1.0
                    else:
                        sim = self.base_pair_comparison_similarity(structures[i], structures[j])
                        similarity_matrix[i][j] = sim
                        similarity_matrix[j][i] = sim
            
            # Convert to distance matrix
            distance_matrix = 1 - similarity_matrix
            condensed_distances = squareform(distance_matrix)
            
            # Perform hierarchical clustering
            linkage_matrix = linkage(condensed_distances, method='ward')
            distance_threshold = 1 - self.cluster_threshold
            clusters = fcluster(linkage_matrix, distance_threshold, criterion='distance')
            
            # Add cluster information to dataframe
            df['cluster'] = clusters
            
            # Add normalized structures for reference
            df['normalized_for_clustering'] = df['structure'].apply(self._normalize_structure_for_display)
            
            # Create cluster summary (now sorted by total_count)
            cluster_summary = self.create_cluster_summary(df)
            
            # Save results
            clustered_file = f"{name}_clustered_structures{pos_suffix}.csv"
            cluster_summary_file = f"{name}_cluster_summary{pos_suffix}.csv"
            
            df.to_csv(clustered_file, index=False)
            cluster_summary.to_csv(cluster_summary_file, index=False)
            
            print(f"✓ Created {clustered_file}")
            print(f"✓ Created {cluster_summary_file}")
            print(f"  Found {len(cluster_summary)} clusters with threshold {self.cluster_threshold}")
            print(f"  Results ordered by total count (most dominant clusters first)")
            
            return clustered_file, cluster_summary_file
            
        except Exception as e:
            print(f"Error in clustering analysis: {e}")
            return None

    def _normalize_structure_for_display(self, structure):
        """Normalize structure for display - shows what clustering algorithm sees"""
        normalized = ""
        for char in structure:
            if char in '([{<':
                normalized += '('
            elif char in ')]}> ':
                normalized += ')'
            else:
                normalized += char
        return normalized

def main():
    parser = argparse.ArgumentParser(description='RNA structure temperature analysis with position-specific options and structure analysis')
    parser.add_argument('--name', '-n', required=True, help='Base name for files')
    parser.add_argument('--min-temp', '-t', type=int, required=True, help='Min temperature level')
    parser.add_argument('--max-temp', '-T', type=int, required=True, help='Max temperature level')
    parser.add_argument('--structure', '-s', required=True, help='Structure file pattern')
    parser.add_argument('--min-frame', '-f', type=int, default=1000, help='Min frame number (default: 1000)')
    parser.add_argument('--positions', '-p', help='Positions to analyze (e.g., "24-29" or "24,25,26,27")')
    parser.add_argument('--position-threshold', type=float, 
                       help='Minimum position accuracy threshold to include frames')
    parser.add_argument('--target-structure', help='Target structure for similarity comparison (file or structure string)')
    parser.add_argument('--continue-on-error', action='store_true', help='Continue if one level fails')
    parser.add_argument('--skip-plots', action='store_true', help='Skip generating plots')
    parser.add_argument('--skip-structure-analysis', action='store_true', help='Skip structure collection and analysis')
    parser.add_argument('--skip-clustering', action='store_true', help='Skip clustering analysis')
    parser.add_argument('--cluster-threshold', type=float, default=0.7, help='Clustering threshold (default: 0.7)')
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
    
    # Load target structure if provided
    target_structure = None
    if args.target_structure:
        if os.path.exists(args.target_structure):
            # Read from file
            with open(args.target_structure, 'r') as f:
                target_structure = ''.join(line.strip() for line in f if line.strip())
        else:
            # Use as direct structure string
            target_structure = args.target_structure
        
        if not args.quiet:
            print(f"Target structure loaded: {target_structure[:50]}{'...' if len(target_structure) > 50 else ''}")
    
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
    
    # Structure analysis setup
    structure_collector = None
    if not args.skip_structure_analysis:
        structure_collector = StructureCollector(target_structure)
    
    # Process temperature levels
    successful_levels = []
    failed_levels = []
    
    for temp_level in range(args.min_temp, args.max_temp + 1):
        success = process_temperature_level(
            temp_level, args.name, args.structure, args.min_frame, script_dir, 
            positions, args.position_threshold, args.quiet,
            collect_structures=not args.skip_structure_analysis,
            structure_collector=structure_collector
        )
        if success:
            successful_levels.append(temp_level)
        else:
            failed_levels.append(temp_level)
            if not args.continue_on_error:
                break
    
    cleanup_pdbs()
    
    print(f"Completed: {len(successful_levels)}/{args.max_temp - args.min_temp + 1} levels")
    if failed_levels:
        print(f"Failed levels: {failed_levels}")
    
    # Structure analysis (NEW FEATURE)
    unique_structures_file = None
    if not args.skip_structure_analysis and successful_levels and structure_collector:
        all_structures, structure_metadata = structure_collector.get_results()
        
        if all_structures:
            result = create_unique_structures_summary(
                all_structures, structure_metadata, args.name, pos_suffix, target_structure
            )
            if result:
                unique_structures_file, all_structures_file = result
    
    # Clustering analysis (NEW FEATURE)
    if not args.skip_clustering and unique_structures_file:
        analyzer = RNAStructureAnalyzer(cluster_threshold=args.cluster_threshold)
        analyzer.perform_clustering_analysis(unique_structures_file, args.name, pos_suffix)
    
    # Update summary file
    if args.summary and successful_levels:
        update_summary_file(args.summary, args.label, successful_levels, 
                          args.threshold, args.median, args.mean, positions_only, pos_suffix)
    
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