#!/usr/bin/env python3

import re
import subprocess
import argparse

# Argument parser for --cutoff and --prefix
parser = argparse.ArgumentParser(description="Extract major clusters above a percentage cutoff.")
parser.add_argument('--cutoff', type=float, default=5, help='Percentage cutoff for major clusters (default: 5)')
parser.add_argument('--prefix', type=str, default='all_thrs15.00A', help='Prefix for cluster files (default: all_thrs15.00A)')
args = parser.parse_args()

input_file = 'log'  # The log file in the current directory
output_file = f"00clusters_over_{args.cutoff}_percent.txt"

pattern = re.compile(r'cluster\s+(\d+)\s+contains\s+\d+\s+\(([\d.]+)%\)\s+structures')

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        match = pattern.search(line)
        if match:
            cluster_num = int(match.group(1))
            percent = float(match.group(2))
            if percent > args.cutoff:
                outfile.write(f"{cluster_num}_ {percent:.2f}\n")
                # Format cluster number as two digits
                cluster_num_str = f"{cluster_num:02d}"
                trafl_file = f"{args.prefix}_clust{cluster_num_str}.trafl"
                # Build the command
                cmd = [
                    "SimRNA_trafl2pdbs",
                    "*_01-000001.pdb",
                    trafl_file,
                    "1",
                    "AA"
                ]
                # Run the command in the shell so that wildcard expands
                subprocess.run(" ".join(cmd), shell=True)