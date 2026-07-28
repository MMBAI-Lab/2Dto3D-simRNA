#!/usr/bin/env python3
"""
plot_sasa_summary.py

Create plots from SASA boxplot summary data collected from multiple directories.

Usage:
    python plot_sasa_summary.py sasa_summary.csv --output sasa_comparison.png
    python plot_sasa_summary.py sasa_summary.csv --metric rel_sasa --group-by nucleotide
"""

import argparse
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.colors as mcolors
from typing import Optional, List

def read_summary_data(csv_path: str) -> pd.DataFrame:
    """Read and validate the SASA summary CSV file"""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Summary file not found: {csv_path}")
    
    df = pd.read_csv(csv_path)
    
    # Validate required columns
    required_cols = ['nucleotide', 'run_id']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    return df

def get_metric_columns(df: pd.DataFrame, metric_type: str) -> List[str]:
    """Get columns for a specific metric type (abs or rel)"""
    suffix = f"_{metric_type}"
    metric_cols = [col for col in df.columns if col.endswith(suffix)]
    return metric_cols

def get_extended_colors(n_colors):
    """
    Generate an extended color palette with many distinct colors
    """
    if n_colors <= 10:
        # Use tab10 for better distinction
        return plt.cm.tab10(np.linspace(0, 1, n_colors))
    elif n_colors <= 20:
        # Use tab20 which has 20 maximally distinct colors
        return plt.cm.tab20(np.linspace(0, 1, n_colors))
    else:
        # For many colors, combine multiple palettes strategically
        base_colors = []
        
        # Start with tab20 (20 very distinct colors)
        base_colors.extend(plt.cm.tab20(np.linspace(0, 1, 20)))
        remaining = n_colors - 20
        
        if remaining > 0:
            # Add Set1 (9 colors, but avoid overlap with tab20)
            set1_colors = plt.cm.Set1(np.linspace(0, 1, min(9, remaining)))
            base_colors.extend(set1_colors)
            remaining -= min(9, remaining)
        
        if remaining > 0:
            # Add Set2 (8 colors)
            set2_colors = plt.cm.Set2(np.linspace(0, 1, min(8, remaining)))
            base_colors.extend(set2_colors)
            remaining -= min(8, remaining)
        
        if remaining > 0:
            # Add Dark2 (8 colors)
            dark2_colors = plt.cm.Dark2(np.linspace(0, 1, min(8, remaining)))
            base_colors.extend(dark2_colors)
            remaining -= min(8, remaining)
        
        if remaining > 0:
            # Add Set3 (12 colors)
            set3_colors = plt.cm.Set3(np.linspace(0, 1, min(12, remaining)))
            base_colors.extend(set3_colors)
            remaining -= min(12, remaining)
        
        if remaining > 0:
            # Use hsv for any remaining colors
            hsv_colors = plt.cm.hsv(np.linspace(0, 1, remaining))
            base_colors.extend(hsv_colors)
        
        return np.array(base_colors[:n_colors])

def create_comparison_boxplot(df: pd.DataFrame, metric_type: str = 'abs', 
                            output_path: Optional[str] = None, 
                            group_by: str = 'run_id',
                            nucleotides: Optional[List[str]] = None) -> str:
    """Create comparison boxplot across runs or nucleotides"""
    
    # Filter data if specific nucleotides requested
    if nucleotides:
        df = df[df['nucleotide'].isin(nucleotides)].copy()
    
    if df.empty:
        raise ValueError("No data remaining after filtering")
    
    # Get metric columns
    metric_cols = get_metric_columns(df, metric_type)
    if not metric_cols:
        raise ValueError(f"No {metric_type} metric columns found")
    
    # Prepare data for plotting
    plot_data = []
    for _, row in df.iterrows():
        base_info = {
            'nucleotide': row['nucleotide'],
            'run_id': row['run_id']
        }
        
        for col in ['median', 'mean', 'q1', 'q3', 'min', 'max']:
            full_col = f"{col}_{metric_type}"
            if full_col in df.columns and pd.notna(row[full_col]):
                plot_data.append({
                    **base_info,
                    'statistic': col,
                    'value': row[full_col]
                })
    
    if not plot_data:
        raise ValueError("No valid data for plotting")
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create the plot
    plt.figure(figsize=(12, 8))
    
    if group_by == 'nucleotide':
        # Group by nucleotide, compare across runs (hue=run_id)
        unique_runs = sorted(df['run_id'].unique())
        colors = get_extended_colors(len(unique_runs))
        color_palette = [mcolors.to_hex(color) for color in colors]
        color_map = dict(zip(unique_runs, color_palette))
        print(f"Coloring by run_id: {color_map}")
        
        ax = sns.boxplot(data=plot_df[plot_df['statistic'].isin(['median', 'q1', 'q3'])], 
                        x='nucleotide', y='value', hue='run_id', palette=color_map)
        plt.title(f'SASA {metric_type.upper()} Comparison Across Runs')
        plt.xlabel('Nucleotide')
    else:
        # Group by run_id, compare across nucleotides (hue=nucleotide)
        unique_nucleotides = sorted(df['nucleotide'].unique())
        
        if len(unique_nucleotides) == 1:
            # Special case: only one nucleotide, so color by run_id instead
            unique_runs = sorted(df['run_id'].unique())
            colors = get_extended_colors(len(unique_runs))
            color_palette = [mcolors.to_hex(color) for color in colors]
            color_map = dict(zip(unique_runs, color_palette))
            print(f"Single nucleotide detected - coloring by run_id: {color_map}")
            
            # Create a custom plot where each run gets its own color
            plot_df_filtered = plot_df[plot_df['statistic'].isin(['median', 'q1', 'q3'])]
            ax = sns.boxplot(data=plot_df_filtered, x='run_id', y='value', palette=color_palette)
            plt.title(f'SASA {metric_type.upper()} Comparison Across Runs')
            plt.xlabel('Run ID')
        else:
            # Multiple nucleotides, color by nucleotide
            colors = get_extended_colors(len(unique_nucleotides))
            color_palette = [mcolors.to_hex(color) for color in colors]
            color_map = dict(zip(unique_nucleotides, color_palette))
            print(f"Coloring by nucleotide: {color_map}")
            
            ax = sns.boxplot(data=plot_df[plot_df['statistic'].isin(['median', 'q1', 'q3'])], 
                            x='run_id', y='value', hue='nucleotide', palette=color_map)
            plt.title(f'SASA {metric_type.upper()} Comparison Across Nucleotides')
            plt.xlabel('Run ID')
        
        plt.xticks(rotation=45, ha='right')
    
    metric_label = 'Relative SASA' if metric_type == 'rel' else 'SASA (Ų)'
    plt.ylabel(metric_label)
    
    if metric_type == 'rel':
        plt.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='Fully exposed')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plot_path = output_path
    else:
        plot_path = f"sasa_{metric_type}_comparison_{group_by}.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    
    plt.close()
    return plot_path

def create_violin_plot(df: pd.DataFrame, metric_type: str = 'abs',
                      output_path: Optional[str] = None,
                      group_by: str = 'run_id',
                      nucleotides: Optional[List[str]] = None) -> str:
    """Create violin plot showing distribution statistics"""
    
    # Filter data if specific nucleotides requested
    if nucleotides:
        df = df[df['nucleotide'].isin(nucleotides)].copy()
    
    # Prepare data - create synthetic distributions from summary statistics
    plot_data = []
    for _, row in df.iterrows():
        median_col = f'median_{metric_type}'
        std_col = f'std_{metric_type}'
        min_col = f'min_{metric_type}'
        max_col = f'max_{metric_type}'
        
        if all(col in df.columns and pd.notna(row[col]) for col in [median_col, std_col]):
            # Create synthetic data points based on statistics
            median_val = row[median_col]
            std_val = row[std_col]
            min_val = row.get(min_col, median_val - 2*std_val) if pd.notna(row.get(min_col)) else median_val - 2*std_val
            max_val = row.get(max_col, median_val + 2*std_val) if pd.notna(row.get(max_col)) else median_val + 2*std_val
            
            # Generate synthetic distribution
            n_points = 100
            synthetic_values = np.random.normal(median_val, std_val, n_points)
            synthetic_values = np.clip(synthetic_values, min_val, max_val)
            
            for val in synthetic_values:
                plot_data.append({
                    'nucleotide': row['nucleotide'],
                    'run_id': row['run_id'],
                    'value': val
                })
    
    if not plot_data:
        raise ValueError("No valid data for violin plot")
    
    plot_df = pd.DataFrame(plot_data)
    
    plt.figure(figsize=(12, 8))
    
    if group_by == 'nucleotide':
        # Group by nucleotide, compare across runs (hue=run_id)
        unique_runs = sorted(df['run_id'].unique())
        colors = get_extended_colors(len(unique_runs))
        color_palette = [mcolors.to_hex(color) for color in colors]
        color_map = dict(zip(unique_runs, color_palette))
        
        ax = sns.violinplot(data=plot_df, x='nucleotide', y='value', hue='run_id', 
                           inner='box', palette=color_map)
        plt.title(f'SASA {metric_type.upper()} Distribution Comparison Across Runs')
        plt.xlabel('Nucleotide')
    else:
        # Group by run_id, compare across nucleotides (hue=nucleotide)
        unique_nucleotides = sorted(df['nucleotide'].unique())
        
        if len(unique_nucleotides) == 1:
            # Special case: only one nucleotide, so color by run_id instead
            unique_runs = sorted(df['run_id'].unique())
            colors = get_extended_colors(len(unique_runs))
            color_palette = [mcolors.to_hex(color) for color in colors]
            
            # Create a custom plot where each run gets its own color
            ax = sns.violinplot(data=plot_df, x='run_id', y='value', palette=color_palette, inner='box')
            plt.title(f'SASA {metric_type.upper()} Distribution Comparison Across Runs')
            plt.xlabel('Run ID')
        else:
            # Multiple nucleotides, color by nucleotide
            colors = get_extended_colors(len(unique_nucleotides))
            color_palette = [mcolors.to_hex(color) for color in colors]
            color_map = dict(zip(unique_nucleotides, color_palette))
            
            ax = sns.violinplot(data=plot_df, x='run_id', y='value', hue='nucleotide', 
                               inner='box', palette=color_map)
            plt.title(f'SASA {metric_type.upper()} Distribution Comparison Across Nucleotides')
            plt.xlabel('Run ID')
        
        plt.xticks(rotation=45, ha='right')
    
    metric_label = 'Relative SASA' if metric_type == 'rel' else 'SASA (Ų)'
    plt.ylabel(metric_label)
    
    if metric_type == 'rel':
        plt.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='Fully exposed')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plot_path = output_path
    else:
        plot_path = f"sasa_{metric_type}_violin_{group_by}.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    
    plt.close()
    return plot_path

def create_heatmap(df: pd.DataFrame, metric: str = 'median_abs',
                  output_path: Optional[str] = None) -> str:
    """Create heatmap showing metric values across runs and nucleotides"""
    
    if metric not in df.columns:
        available_metrics = [col for col in df.columns if any(x in col for x in ['median', 'mean', 'std', 'min', 'max'])]
        raise ValueError(f"Metric '{metric}' not found. Available: {available_metrics}")
    
    # Pivot data for heatmap
    pivot_df = df.pivot(index='run_id', columns='nucleotide', values=metric)
    
    plt.figure(figsize=(10, 8))
    
    # Create heatmap
    ax = sns.heatmap(pivot_df, annot=True, fmt='.2f', cmap='viridis', 
                     cbar_kws={'label': metric.replace('_', ' ').title()})
    
    plt.title(f'SASA {metric.replace("_", " ").title()} Heatmap')
    plt.xlabel('Nucleotide')
    plt.ylabel('Run ID')
    plt.tight_layout()
    
    # Save plot
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plot_path = output_path
    else:
        plot_path = f"sasa_{metric}_heatmap.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    
    plt.close()
    return plot_path

def create_summary_stats_table(df: pd.DataFrame, output_path: Optional[str] = None) -> str:
    """Create a summary statistics table"""
    
    # Calculate summary statistics
    summary_data = []
    
    for nucleotide in df['nucleotide'].unique():
        nuc_data = df[df['nucleotide'] == nucleotide]
        
        for metric_type in ['abs', 'rel']:
            median_col = f'median_{metric_type}'
            if median_col in df.columns:
                median_vals = nuc_data[median_col].dropna()
                if not median_vals.empty:
                    summary_data.append({
                        'Nucleotide': nucleotide,
                        'Metric': f'Median {metric_type.upper()}',
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
        table_path = "sasa_summary_stats.csv"
        summary_df.to_csv(table_path, index=False, float_format='%.4f')
    
    return table_path

def main():
    parser = argparse.ArgumentParser(description='Plot SASA summary data from multiple directories')
    parser.add_argument('summary_csv', help='Path to SASA summary CSV file')
    parser.add_argument('--output', '-o', help='Output plot path (default: auto-generate)')
    parser.add_argument('--metric', choices=['abs', 'rel'], default='rel', 
                       help='Metric type to plot: abs (absolute SASA) or rel (relative SASA)')
    parser.add_argument('--plot-type', choices=['boxplot', 'violin', 'heatmap'], default='boxplot',
                       help='Type of plot to create')
    parser.add_argument('--group-by', choices=['run_id', 'nucleotide'], default='run_id',
                       help='How to group the data for comparison')
    parser.add_argument('--nucleotides', nargs='+', 
                       help='Specific nucleotides to include (e.g., C:1 A:10)')
    parser.add_argument('--heatmap-metric', default='median_abs',
                       help='Specific metric for heatmap (e.g., median_abs, mean_rel)')
    parser.add_argument('--create-all', action='store_true',
                       help='Create all plot types and summary table')
    parser.add_argument('--summary-table', action='store_true',
                       help='Create summary statistics table')
    
    args = parser.parse_args()
    
    try:
        # Read data
        print(f"Reading summary data from {args.summary_csv}...")
        df = read_summary_data(args.summary_csv)
        print(f"✓ Loaded {len(df)} rows with {len(df['run_id'].unique())} runs and {len(df['nucleotide'].unique())} nucleotides")
        
        # Check available metrics
        abs_metrics = get_metric_columns(df, 'abs')
        rel_metrics = get_metric_columns(df, 'rel')
        print(f"Available metrics: {len(abs_metrics)} absolute, {len(rel_metrics)} relative")
        
        if args.metric == 'rel' and not rel_metrics:
            print("Warning: Relative SASA metrics not found, switching to absolute")
            args.metric = 'abs'
        
        created_files = []
        
        if args.create_all or args.plot_type == 'boxplot':
            print("Creating boxplot...")
            try:
                plot_path = create_comparison_boxplot(df, args.metric, args.output, 
                                                    args.group_by, args.nucleotides)
                created_files.append(plot_path)
                print(f"✓ Created: {plot_path}")
            except Exception as e:
                print(f"Error creating boxplot: {e}")
        
        if args.create_all or args.plot_type == 'violin':
            print("Creating violin plot...")
            try:
                plot_path = create_violin_plot(df, args.metric, args.output, 
                                             args.group_by, args.nucleotides)
                created_files.append(plot_path)
                print(f"✓ Created: {plot_path}")
            except Exception as e:
                print(f"Error creating violin plot: {e}")
        
        if args.create_all or args.plot_type == 'heatmap':
            print("Creating heatmap...")
            try:
                plot_path = create_heatmap(df, args.heatmap_metric, args.output)
                created_files.append(plot_path)
                print(f"✓ Created: {plot_path}")
            except Exception as e:
                print(f"Error creating heatmap: {e}")
        
        if args.create_all or args.summary_table:
            print("Creating summary table...")
            try:
                table_path = create_summary_stats_table(df, args.output)
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