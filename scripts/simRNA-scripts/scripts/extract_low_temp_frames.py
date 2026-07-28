#!/usr/bin/env python3
"""
Script to extract frames from TRAFL files based on temperature levels from log files.
Optimized version with caching, indexing, and multiprocessing.
"""

import argparse
import re
import os
import sys
import glob
import time
from typing import Dict, Set, Tuple, List
from multiprocessing import Pool, cpu_count
from functools import partial

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    print("Note: Install tqdm for progress bars: pip install tqdm")

class TraflIndex:
    """Index for fast frame lookup in TRAFL files"""
    
    def __init__(self, trafl_file: str):
        self.trafl_file = trafl_file
        self.frame_positions = {}  # frame_number -> (line_start, line_end)
        self.file_lines = None
        self._build_index()
    
    def _build_index(self):
        """Build index of frame positions"""
        try:
            with open(self.trafl_file, 'r') as f:
                self.file_lines = f.readlines()
        except FileNotFoundError:
            return
        
        for i, line in enumerate(self.file_lines):
            line_stripped = line.strip()
            if line_stripped and not line_stripped.startswith('#'):
                parts = line_stripped.split()
                if parts and parts[0].isdigit():
                    frame_number = int(parts[0])
                    # Each frame is 2 lines (current + next)
                    self.frame_positions[frame_number] = (i, min(i + 2, len(self.file_lines)))
    
    def get_frame(self, frame_number: int) -> str:
        """Get frame content by number"""
        if frame_number not in self.frame_positions:
            return ""
        
        start, end = self.frame_positions[frame_number]
        return ''.join(self.file_lines[start:end])

def parse_log_file_optimized(log_file: str, min_temp: float, max_temp: float, min_frame: int) -> Dict[int, Set[int]]:
    """
    Optimized log file parsing with better regex compilation.
    """
    replica_frames = {}
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"Error: Log file '{log_file}' not found")
        sys.exit(1)
    
    print("Parsing log file...")
    
    # Compile regex patterns once
    temp_pattern = re.compile(r'temp\.\s*level:\s*(\d+),\s*replica:\s*(\d+)')
    write_pattern = re.compile(r'Write number:\s*(\d+)')
    
    # Split by separator and process in chunks
    sections = content.split('==================================')
    
    if HAS_TQDM:
        sections_iter = tqdm(range(len(sections) - 1), desc="Processing log sections")
    else:
        sections_iter = range(len(sections) - 1)
        print(f"Processing {len(sections) - 1} log sections...")
    
    for i in sections_iter:
        section = sections[i] + sections[i + 1]
        
        # Use compiled patterns
        temp_match = temp_pattern.search(section)
        if temp_match:
            temp_level = int(temp_match.group(1))
            replica = int(temp_match.group(2))
            
            # Check if temperature is in desired range
            if min_temp <= temp_level <= max_temp:
                write_match = write_pattern.search(section)
                if write_match:
                    frame_number = int(write_match.group(1))
                    
                    # Check if frame number meets minimum requirement
                    if frame_number >= min_frame:
                        if replica not in replica_frames:
                            replica_frames[replica] = set()
                        replica_frames[replica].add(frame_number)
    
    return replica_frames

def extract_frames_from_replica(args: Tuple[int, Set[int], str]) -> List[str]:
    """
    Extract frames from a single replica (for multiprocessing).
    """
    replica, frames, base_name = args
    trafl_file = f"{base_name}_{replica:02d}.trafl"
    
    # Build index for this file
    index = TraflIndex(trafl_file)
    if not index.file_lines:
        return []
    
    extracted = []
    for frame_number in sorted(frames):
        frame_content = index.get_frame(frame_number)
        if frame_content:
            extracted.append(frame_content)
    
    return extracted

def main():
    parser = argparse.ArgumentParser(description='Extract low temperature frames from TRAFL files')
    parser.add_argument('log_file', help='Path to the C*.log file')
    parser.add_argument('--min-temp', type=float, default=0, help='Minimum temperature level (default: 0)')
    parser.add_argument('--max-temp', type=float, default=10, help='Maximum temperature level (default: 10)')
    parser.add_argument('--min-frame', type=int, default=0, help='Minimum frame number to consider (default: 0)')
    parser.add_argument('--base-name', required=True, help='Base name for TRAFL files (e.g., "name" for name_01.trafl)')
    parser.add_argument('--output', default='low_temp.trafl', help='Output compilation file (default: low_temp.trafl)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    parser.add_argument('--processes', type=int, default=None, help='Number of processes (default: auto)')
    parser.add_argument('--single-threaded', action='store_true', help='Disable multiprocessing')
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    # Expand glob pattern if needed
    if '*' in args.log_file:
        log_files = glob.glob(args.log_file)
        if not log_files:
            print(f"No files match pattern: {args.log_file}")
            sys.exit(1)
        args.log_file = log_files[0]  # Use first match
        if args.verbose:
            print(f"Using log file: {args.log_file}")
    
    print(f"Extracting frames from temperature levels {args.min_temp}-{args.max_temp}, min frame {args.min_frame}")
    
    # Parse the log file to get replica-frame pairs
    replica_frames = parse_log_file_optimized(args.log_file, args.min_temp, args.max_temp, args.min_frame)
    
    if not replica_frames:
        print("No frames found matching the criteria")
        sys.exit(0)
    
    # Count total unique frames
    total_frames = sum(len(frames) for frames in replica_frames.values())
    print(f"Found {total_frames} unique frames across {len(replica_frames)} replicas")
    
    if args.verbose:
        for replica, frames in replica_frames.items():
            print(f"  Replica {replica}: {len(frames)} unique frames")
    
    # Determine number of processes
    if args.single_threaded:
        num_processes = 1
    else:
        num_processes = args.processes or min(cpu_count(), len(replica_frames))
    
    print(f"Using {num_processes} process(es) for extraction")
    
    # Prepare arguments for multiprocessing
    extraction_args = [(replica, frames, args.base_name) for replica, frames in replica_frames.items()]
    
    # Extract frames using multiprocessing
    print("Extracting frames...")
    if num_processes == 1:
        # Single-threaded execution
        if HAS_TQDM:
            all_extracted = []
            for arg in tqdm(extraction_args, desc="Processing replicas"):
                all_extracted.extend(extract_frames_from_replica(arg))
        else:
            all_extracted = []
            for i, arg in enumerate(extraction_args):
                all_extracted.extend(extract_frames_from_replica(arg))
                print(f"Processed replica {i+1}/{len(extraction_args)}")
    else:
        # Multi-threaded execution
        with Pool(processes=num_processes) as pool:
            if HAS_TQDM:
                results = list(tqdm(
                    pool.imap(extract_frames_from_replica, extraction_args),
                    total=len(extraction_args),
                    desc="Processing replicas"
                ))
            else:
                results = pool.map(extract_frames_from_replica, extraction_args)
                print(f"Processed {len(results)} replicas")
        
        # Flatten results
        all_extracted = []
        for result in results:
            all_extracted.extend(result)
    
    # Write compilation file
    if all_extracted:
        print(f"Writing {len(all_extracted)} frames to output file...")
        try:
            with open(args.output, 'w') as f:
                if HAS_TQDM:
                    for frame in tqdm(all_extracted, desc="Writing frames"):
                        f.write(frame)
                else:
                    for i, frame in enumerate(all_extracted):
                        f.write(frame)
                        if i % 5000 == 0 and i > 0:
                            print(f"Written {i}/{len(all_extracted)} frames")
            
            elapsed_time = time.time() - start_time
            print(f"Successfully wrote {len(all_extracted)} frames to {args.output}")
            print(f"Total time: {elapsed_time:.2f} seconds")
            
            if elapsed_time > 0:
                rate = len(all_extracted) / elapsed_time
                print(f"Processing rate: {rate:.1f} frames/second")
                
        except IOError as e:
            print(f"Error writing output file: {e}")
            sys.exit(1)
    else:
        print("No frames were successfully extracted")

if __name__ == "__main__":
    main()