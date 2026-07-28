#!/usr/bin/env python3
"""
plot_summary.py

Create plots from boxplot summary data collected from multiple directories.
Works with both SASA and distance analysis summary files.

Usage:
    python plot_summary.py summary.csv --output comparison.png
    python plot_summary.py sasa_summary.csv --metric rel_sasa --group-by nucleotide
    python plot_summary.py distance_summary.csv --metric distance --group-by run_id
"""

import argparse
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.colors as mcolors
from typing import Optional, List, Dict, Tuple

def detect_data_type(df: pd.DataFrame) -> str:
    """Detect whether this is SASA or distance data based on columns"""
    sasa_indicators = ['sasa_A2', 'rel_sasa', 'nucleotide', 'median_abs', 'median_rel']
    distance_indicators = ['distance_A', 'nt1', 'nt2', 'median', 'mean', 'std']
    
    sasa_score = sum(1 for col in df.columns if any(ind in col.lower() for ind in sasa_indicators))
    distance_score = sum(1 for col in df.columns if any(ind in col.lower() for ind in distance_indicators))
    
    if sasa_score > distance_score:
        return 'sasa'
    elif 'nt1' in df.columns and 'nt2' in df.columns:
        return 'distance'
    elif 'nucleotide' in df.columns:
        return 'sasa'
    else:
        return 'unknown'

def get_data_info(df: pd.DataFrame) -> Dict:
    """Get information about the data structure"""
    data_type = detect_data_type(df)
    info = {'type': data_type, 'metrics': [], 'group_columns': []}
    
    if data_type == 'sasa':
        # SASA data structure
        abs_metrics = [col for col in df.columns if col.endswith('_abs')]
        rel_metrics = [col for col in df.columns if col.endswith('_rel')]
        info['metrics'] = {
            'abs': abs_metrics,
            'rel': rel_metrics
        }
        info['group_columns'] = ['nucleotide', 'run_id']
        info['primary_id'] = 'nucleotide'
        info['units'] = {'abs': 'SASA (Ų)', 'rel': 'Relative SASA'}
        
    elif data_type == 'distance':
        # Distance data structure
        basic_metrics = [col for col in df.columns if col in ['median', 'mean', 'std', 'min', 'max', 'q1', 'q3']]
        info['metrics'] = {'distance': basic_metrics}
        info['group_columns'] = ['nt1', 'nt2', 'run_id']
        info['primary_id'] = 'pair'  # We'll create this from nt1+nt2
        info['units'] = {'distance': 'Distance (Å)'}
        
        # Create pair identifier
        if 'nt1' in df.columns and 'nt2' in df.columns:
            df['pair'] = df['nt1'].astype(str) + ' — ' + df['nt2'].astype(str)
    
    return info

def read_summary_data(csv_path: str) -> Tuple[pd.DataFrame, Dict]:
    """Read and analyze the summary CSV file"""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Summary file not found: {csv_path}")
    
    df = pd.read_csv(csv_path)
    
    # Validate required columns
    if 'run_id' not in df.columns:
        raise ValueError("Missing required column: run_id")
    
    info = get_data_info(df)
    print(f"Detected data type: {info['type']}")
    print(f"Available metrics: {info['metrics']}")
    
    return df, info

def get_extended_colors(n_colors):
    """Generate an extended color palette with many distinct colors"""
    if n_colors <= 10:
        return plt.cm.tab10(np.linspace(0, 1, n_colors))
    elif n_colors <= 20:
        return plt.cm.tab20(np.linspace(0, 1, n_colors))
    else:
        base_colors = []
        base_colors.extend(plt.cm.tab20(np.linspace(0, 1, 20)))
        remaining = n_colors - 20
        
        if remaining > 0:
            set1_colors = plt.cm.Set1(np.linspace(0, 1, min(9, remaining)))
            base_colors.extend(set1_colors)
            remaining -= min(9, remaining)
        
        if remaining > 0:
            set2_colors = plt.cm.Set2(np.linspace(0, 1, min(8, remaining)))
            base_colors.extend(set2_colors)
            remaining -= min(8, remaining)
        
        if remaining > 0:
            hsv_colors = plt.cm.hsv(np.linspace(0, 1, remaining))
            base_colors.extend(hsv_colors)
        
        return np.array(base_colors[:n_colors])

def create_comparison_boxplot(df: pd.DataFrame, info: Dict, metric: str = None,
                            output_path: Optional[str] = None, 
                            group_by: str = 'run_id',
                            filter_items: Optional[List[str]] = None) -> str:
    """Create comparison boxplot across runs or items"""
    
    # Auto-select metric if not provided
    if metric is None:
        if info['type'] == 'sasa':
            metric = 'rel' if info['metrics']['rel'] else 'abs'
        elif info['type'] == 'distance':
            metric = 'distance'
    
    # Get the appropriate metric column
    if info['type'] == 'sasa':
        if metric in ['rel', 'relative']:
            metric_col = 'median_rel'
            ylabel = info['units']['rel']
        else:
            metric_col = 'median_abs' 
            ylabel = info['units']['abs']
    elif info['type'] == 'distance':
        metric_col = 'median'
        ylabel = info['units']['distance']
    else:
        raise ValueError(f"Unknown data type: {info['type']}")
    
    if metric_col not in df.columns:
        available = [col for col in df.columns if 'median' in col]
        raise ValueError(f"Metric column '{metric_col}' not found. Available: {available}")
    
    # Filter data if specific items requested
    if filter_items:
        if info['type'] == 'sasa':
            df = df[df['nucleotide'].isin(filter_items)].copy()
        elif info['type'] == 'distance':
            df = df[df['pair'].isin(filter_items)].copy()
    
    if df.empty:
        raise ValueError("No data remaining after filtering")
    
    # Prepare data for plotting - create synthetic distributions from summary stats
    plot_data = []
    for _, row in df.iterrows():
        if pd.notna(row[metric_col]):
            # Get related statistics
            if info['type'] == 'sasa':
                suffix = '_rel' if 'rel' in metric_col else '_abs'
                std_col = f'std{suffix}'
                min_col = f'min{suffix}'
                max_col = f'max{suffix}'
                q1_col = f'q1{suffix}'
                q3_col = f'q3{suffix}'
            else:
                std_col = 'std'
                min_col = 'min'
                max_col = 'max'
                q1_col = 'q1'
                q3_col = 'q3'
            
            median_val = row[metric_col]
            std_val = row.get(std_col, 0) if pd.notna(row.get(std_col)) else 0
            q1_val = row.get(q1_col, median_val - std_val) if pd.notna(row.get(q1_col)) else median_val - std_val
            q3_val = row.get(q3_col, median_val + std_val) if pd.notna(row.get(q3_col)) else median_val + std_val
            
            # Add the actual statistics as data points
            if pd.notna(q1_val):
                plot_data.append({
                    info['primary_id']: row.get(info['primary_id'], row.get('nucleotide', row.get('pair', 'unknown'))),
                    'run_id': row['run_id'],
                    'value': q1_val,
                    'statistic': 'q1'
                })
            
            plot_data.append({
                info['primary_id']: row.get(info['primary_id'], row.get('nucleotide', row.get('pair', 'unknown'))),
                'run_id': row['run_id'],
                'value': median_val,
                'statistic': 'median'
            })
            
            if pd.notna(q3_val):
                plot_data.append({
                    info['primary_id']: row.get(info['primary_id'], row.get('nucleotide', row.get('pair', 'unknown'))),
                    'run_id': row['run_id'],
                    'value': q3_val,
                    'statistic': 'q3'
                })
    
    if not plot_data:
        raise ValueError("No valid data for plotting")
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    primary_col = info['primary_id']
    if group_by == primary_col or group_by == 'nucleotide' or group_by == 'pair':
        # Group by primary identifier (nucleotide/pair), compare across runs
        unique_runs = sorted(df['run_id'].unique())
        colors = get_extended_colors(len(unique_runs))
        color_palette = [mcolors.to_hex(color) for color in colors]
        color_map = dict(zip(unique_runs, color_palette))
        print(f"Coloring by run_id: {color_map}")
        
        ax = sns.boxplot(data=plot_df[plot_df['statistic'].isin(['median', 'q1', 'q3'])], 
                        x=primary_col, y='value', hue='run_id', palette=color_map)
        plt.title(f'{info["type"].upper()} {metric} Comparison Across Runs')
        plt.xlabel(primary_col.title())
    else:
        # Group by run_id, compare across primary identifiers
        unique_items = sorted(df[primary_col].unique()) if primary_col in df.columns else []
        
        if len(unique_items) <= 1:
            # Special case: only one item, so color by run_id instead
            unique_runs = sorted(df['run_id'].unique())
            colors = get_extended_colors(len(unique_runs))
            color_palette = [mcolors.to_hex(color) for color in colors]
            
            ax = sns.boxplot(data=plot_df[plot_df['statistic'].isin(['median', 'q1', 'q3'])], 
                            x='run_id', y='value', palette=color_palette)
            plt.title(f'{info["type"].upper()} {metric} Comparison Across Runs')
            plt.xlabel('Run ID')
        else:
            # Multiple items, color by primary identifier
            colors = get_extended_colors(len(unique_items))
            color_palette = [mcolors.to_hex(color) for color in colors]
            color_map = dict(zip(unique_items, color_palette))
            print(f"Coloring by {primary_col}: {color_map}")
            
            ax = sns.boxplot(data=plot_df[plot_df['statistic'].isin(['median', 'q1', 'q3'])], 
                            x='run_id', y='value', hue=primary_col, palette=color_map)
            plt.title(f'{info["type"].upper()} {metric} Comparison Across {primary_col.title()}s')
            plt.xlabel('Run ID')
        
        plt.xticks(rotation=45, ha='right')
    
    plt.ylabel(ylabel)
    
    # Add reference lines for specific metrics
    if info['type'] == 'sasa' and metric in ['rel', 'relative']:
        plt.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='Fully exposed')
        plt.legend()
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plot_path = output_path
    else:
        plot_path = f"{info['type']}_{metric}_comparison_{group_by}.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    
    plt.close()
    return plot_path

def create_heatmap(df: pd.DataFrame, info: Dict, metric: str = None,
                  output_path: Optional[str] = None) -> str:
    """Create heatmap showing metric values across runs and items"""
    
    # Auto-select metric if not provided
    if metric is None:
        if info['type'] == 'sasa':
            metric_col = 'median_rel' if 'median_rel' in df.columns else 'median_abs'
        else:
            metric_col = 'median'
    else:
        if info['type'] == 'sasa':
            if metric in ['rel', 'relative']:
                metric_col = 'median_rel'
            else:
                metric_col = 'median_abs'
        else:
            metric_col = 'median'
    
    if metric_col not in df.columns:
        available_metrics = [col for col in df.columns if 'median' in col]
        raise ValueError(f"Metric '{metric_col}' not found. Available: {available_metrics}")
    
    primary_col = info['primary_id']
    
    # Pivot data for heatmap
    pivot_df = df.pivot(index='run_id', columns=primary_col, values=metric_col)
    
    plt.figure(figsize=(10, 8))
    
    # Create heatmap
    ax = sns.heatmap(pivot_df, annot=True, fmt='.2f', cmap='viridis', 
                     cbar_kws={'label': metric_col.replace('_', ' ').title()})
    
    plt.title(f'{info["type"].upper()} {metric_col.replace("_", " ").title()} Heatmap')
    plt.xlabel(primary_col.title())
    plt.ylabel('Run ID')
    plt.tight_layout()
    
    # Save plot
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plot_path = output_path
    else:
        plot_path = f"{info['type']}_{metric_col}_heatmap.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    
    plt.close()
    return plot_path

def create_summary_stats_table(df: pd.DataFrame, info: Dict, 
                              output_path: Optional[str] = None) -> str:
    """Create a summary statistics table"""
    
    summary_data = []
    primary_col = info['primary_id']
    
    for item in df[primary_col].unique():
        item_data = df[df[primary_col] == item]
        
        if info['type'] == 'sasa':
            for metric_type in ['abs', 'rel']:
                median_col = f'median_{metric_type}'
                if median_col in df.columns:
                    median_vals = item_data[median_col].dropna()
                    if not median_vals.empty:
                        summary_data.append({
                            primary_col.title(): item,
                            'Metric': f'Median {metric_type.upper()}',
                            'Mean': median_vals.mean(),
                            'Std': median_vals.std(),
                            'Min': median_vals.min(),
                            'Max': median_vals.max(),
                            'N_Runs': len(median_vals)
                        })
        else:
            # Distance data
            median_col = 'median'
            if median_col in df.columns:
                median_vals = item_data[median_col].dropna()
                if not median_vals.empty:
                    summary_data.append({
                        primary_col.title(): item,
                        'Metric': 'Distance',
                        'Mean': median_vals.mean(),
                        'Std': median_vals.std(),
                        'Min': median_vals.min(),
                        'Max': median_vals.max(),
                        'N_Runs': len(median_vals)
                    })
    
    if not summary_data:
        raise ValueError("No data for summary table")
    
    summary_df = pd.DataFrame(summary_data)
    
    # Save as CSV
    if output_path:
        summary_df.to_csv(output_path, index=False, float_format='%.4f')
        table_path = output_path
    else:
        table_path = f"{info['type']}_summary_stats.csv"
        summary_df.to_csv(table_path, index=False, float_format='%.4f')
    
    return table_path

def main():
    parser = argparse.ArgumentParser(description='Plot summary data from SASA or distance analyses')
    parser.add_argument('summary_csv', help='Path to summary CSV file')
    parser.add_argument('--output', '-o', help='Output plot path (default: auto-generate)')
    parser.add_argument('--metric', choices=['abs', 'rel', 'relative', 'distance'], 
                       help='Metric type to plot (auto-detected if not specified)')
    parser.add_argument('--plot-type', choices=['boxplot', 'heatmap'], default='boxplot',
                       help='Type of plot to create')
    parser.add_argument('--group-by', choices=['run_id', 'nucleotide', 'pair'], default='run_id',
                       help='How to group the data for comparison')
    parser.add_argument('--filter-items', nargs='+', 
                       help='Specific items to include (nucleotides or pairs)')
    parser.add_argument('--create-all', action='store_true',
                       help='Create all plot types and summary table')
    parser.add_argument('--summary-table', action='store_true',
                       help='Create summary statistics table')
    
    args = parser.parse_args()
    
    try:
        # Read data
        print(f"Reading summary data from {args.summary_csv}...")
        df, info = read_summary_data(args.summary_csv)
        print(f"✓ Loaded {len(df)} rows with {len(df['run_id'].unique())} runs")
        
        if info['type'] == 'sasa':
            print(f"  - {len(df['nucleotide'].unique())} nucleotides")
        elif info['type'] == 'distance':
            print(f"  - {len(df['pair'].unique())} distance pairs")
        
        created_files = []
        
        if args.create_all or args.plot_type == 'boxplot':
            print("Creating boxplot...")
            try:
                plot_path = create_comparison_boxplot(df, info, args.metric, args.output, 
                                                    args.group_by, args.filter_items)
                created_files.append(plot_path)
                print(f"✓ Created: {plot_path}")
            except Exception as e:
                print(f"Error creating boxplot: {e}")
        
        if args.create_all or args.plot_type == 'heatmap':
            print("Creating heatmap...")
            try:
                plot_path = create_heatmap(df, info, args.metric, args.output)
                created_files.append(plot_path)
                print(f"✓ Created: {plot_path}")
            except Exception as e:
                print(f"Error creating heatmap: {e}")
        
        if args.create_all or args.summary_table:
            print("Creating summary table...")
            try:
                table_path = create_summary_stats_table(df, info, args.output)
                created_files.append(table_path)
                print(f"✓ Created: {table_path}")
            except Exception as e:
                print(f"Error creating summary table: {e}")
        
        if created_files:
            print(f"\n✓ Successfully created {len(created_files)} files:")
            for f in created_files:
                print(f"  - {f}")
        else:
            print("No files created")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()