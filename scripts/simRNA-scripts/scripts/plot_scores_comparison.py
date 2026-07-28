#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import argparse

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

def plot_temperature_analysis(csv_file, metric_type=None):
    """
    Create three separate plots: one for pairs data, one for dots data, and one for dots_only data
    with a shared legend and custom metric labeling
    """
    # Read the CSV file
    df = pd.read_csv(csv_file)
    
    # Extract base labels and scoring type
    df['base_label'] = df['label'].str.replace('_pairs|_dots|_dots_only', '', regex=True)
    df['scoring_type'] = df['label'].str.extract(r'_(pairs|dots|dots_only)$')[0]
    
    # Get temperature level columns (numeric columns)
    temp_columns = [col for col in df.columns if col.isdigit()]
    temp_levels = sorted([int(col) for col in temp_columns])
    
    # Get unique base labels and sort them
    base_labels = sorted(df['base_label'].unique())
    
    # Generate extended color palette
    colors = get_extended_colors(len(base_labels))
    
    print(f"Plotting {len(base_labels)} labels with distinct colors")
    print(f"Base labels: {base_labels}")
    
    # Check which scoring types are available
    available_scoring_types = df['scoring_type'].dropna().unique()
    print(f"Available scoring types: {available_scoring_types}")
    
    # Determine number of subplots needed
    n_plots = len(available_scoring_types)
    
    # Create subplots with more space for the shared legend
    fig, axes = plt.subplots(1, n_plots, figsize=(8 * n_plots + 3, 10))
    
    # If only one plot, make axes a list for consistency
    if n_plots == 1:
        axes = [axes]
    
    # Plot each scoring type
    plot_order = ['pairs', 'dots', 'dots_only'] 
    plot_titles = {
        'pairs': 'Base Pairs Only', 
        'dots': 'Including Dots', 
        'dots_only': 'Dots Only'
    }
    
    # Store legend handles and labels from the first plot
    legend_handles = None
    legend_labels = None
    
    plot_idx = 0
    for scoring_type in plot_order:
        if scoring_type in available_scoring_types:
            scoring_data = df[df['scoring_type'] == scoring_type]
            handles, labels = plot_scoring_type(scoring_data, base_labels, temp_levels, colors, 
                                              axes[plot_idx], plot_titles[scoring_type], 
                                              show_legend=(plot_idx == 0), metric_type=metric_type)
            
            # Store legend info from the first plot
            if plot_idx == 0:
                legend_handles = handles
                legend_labels = labels
            
            plot_idx += 1
    
    # Create a single shared legend
    if legend_handles and legend_labels:
        # Determine the best position and formatting for the legend
        if len(base_labels) > 15:
            # For many labels, use smaller font and multiple columns
            ncols = min(3, (len(base_labels) + 9) // 10)  # Max 3 columns
            legend = fig.legend(legend_handles, legend_labels, 
                              bbox_to_anchor=(1.02, 0.5), loc='center left',
                              fontsize=8, ncol=ncols, frameon=True, 
                              fancybox=True, shadow=True)
        elif len(base_labels) > 8:
            # For moderate number of labels, use 2 columns
            legend = fig.legend(legend_handles, legend_labels,
                              bbox_to_anchor=(1.02, 0.5), loc='center left',
                              fontsize=9, ncol=2, frameon=True)
        else:
            # For few labels, single column is fine
            legend = fig.legend(legend_handles, legend_labels,
                              bbox_to_anchor=(1.02, 0.5), loc='center left',
                              fontsize=10, frameon=True)
    
    # Adjust layout to prevent overlap and accommodate the legend
    plt.tight_layout()
    plt.subplots_adjust(right=0.85 if len(base_labels) <= 8 else 0.75)
    
    # Save the plots with metric type in filename if specified
    if metric_type:
        output_file = csv_file.replace('.csv', f'_temperature_analysis_all_three_{metric_type}.png')
    else:
        output_file = csv_file.replace('.csv', '_temperature_analysis_all_three.png')
    
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    print(f"Plots saved as: {output_file}")
    
    # Show the plots
    plt.show()
    
    # Print summary statistics for all scoring types
    print_summary_statistics(df, base_labels, temp_levels, metric_type)

def plot_scoring_type(data, base_labels, temp_levels, colors, ax, title, show_legend=False, metric_type=None):
    """
    Plot data for a specific scoring type (pairs, dots, or dots_only)
    Returns legend handles and labels if show_legend is True
    """
    legend_handles = []
    legend_labels = []
    
    # Plot each base label with its assigned color
    present_labels = [bl for bl in base_labels if bl in data['base_label'].values]
    present_colors = get_extended_colors(len(present_labels))
    for i, base_label in enumerate(present_labels):
        label_data = data[data['base_label'] == base_label]
  
        
        # Calculate statistics for each temperature level
        stats = []
        for temp in temp_levels:
            values = label_data[str(temp)].dropna()

            stats.append({
                'temp': temp,
                'mean': values.mean(),
                'min': values.min(),
                'max': values.max(),
                'count': len(values)
            })

        
        # Convert to arrays for plotting
        temps = [s['temp'] for s in stats]
        means = [s['mean'] for s in stats]
        mins = [s['min'] for s in stats]
        maxs = [s['max'] for s in stats]
        
        # Use the specific color for this label index
        color = present_colors[i]
        
        # Debug: print color assignment
        print(f"Label '{base_label}' (index {i}) assigned color: {color}")
        
        # Plot mean line
        line = ax.plot(temps, means, 'o-', color=color, linewidth=2, 
                      label=f'{base_label}', markersize=6)[0]
        
        # Store legend info from the first subplot
        if show_legend:
            legend_handles.append(line)
            legend_labels.append(base_label)
        
        # Plot min/max as filled area with consistent alpha
        # Use a slightly lower alpha to prevent accumulation effects
        ax.fill_between(temps, mins, maxs, alpha=0.15, color=color, 
                       linewidth=0, edgecolor='none')
        
        # Add min/max as dotted lines with consistent styling
        ax.plot(temps, mins, '--', color=color, alpha=0.6, linewidth=1)
        ax.plot(temps, maxs, '--', color=color, alpha=0.6, linewidth=1)
    
    # Customize the plot
    ax.set_xlabel('Temperature Level', fontsize=12)
    ax.set_title(f'Temperature Analysis: {title}\n(Shaded areas show min-max range)', 
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Set x-axis to show all temperature levels
    ax.set_xticks(temp_levels)
    ax.set_xlim(min(temp_levels) - 0.5, max(temp_levels) + 0.5)
    
    # Determine y-axis label and limits based on metric type or data
    ylabel = determine_ylabel(data, base_labels, temp_levels, metric_type)
    ax.set_ylabel(ylabel, fontsize=12)
    
    # Set y-axis limits based on data type
    all_values = []
    for i, base_label in enumerate(base_labels):
        label_data = data[data['base_label'] == base_label]
        for temp in temp_levels:
            values = label_data[str(temp)].dropna()
            all_values.extend(values.tolist())
    
    if all_values:
        # Check if values look like percentages (0-100) or accuracy (0-1)
        max_val = max(all_values)
        if max_val > 1:
            # Likely percentage values (0-100)
            ax.set_ylim(0, 100)
        else:
            # Likely accuracy values (0-1)
            ax.set_ylim(0, 1)
    
    return legend_handles, legend_labels

def determine_ylabel(data, base_labels, temp_levels, metric_type):
    """
    Determine the appropriate y-axis label based on metric type or data inspection
    """
    if metric_type:
        # Use provided metric type
        metric_lower = metric_type.lower()
        if metric_lower == 'mean':
            return 'Mean Accuracy'
        elif metric_lower == 'median':
            return 'Median Accuracy'
        elif metric_lower == 'percentage':
            return 'Percentage Above Threshold (%)'
        elif metric_lower == 'accuracy':
            return 'Accuracy'
        else:
            # Use the provided metric type as-is with proper capitalization
            return f'{metric_type.title()} Accuracy'
    else:
        # Auto-detect based on data values
        sample_values = []
        for i, base_label in enumerate(base_labels[:3]):  # Check first few labels
            label_data = data[data['base_label'] == base_label]
            for temp in temp_levels[:3]:  # Check first few temperature levels
                values = label_data[str(temp)].dropna()
                sample_values.extend(values.tolist())
        
        if sample_values:
            max_val = max(sample_values)
            if max_val > 1:
                return 'Percentage Above Threshold (%)'
            else:
                return 'Accuracy'
        else:
            return 'Value'

def print_summary_statistics(df, base_labels, temp_levels, metric_type=None):
    """
    Print summary statistics for all three scoring types
    """
    print("\nSummary Statistics:")
    print("=" * 80)
    
    if metric_type:
        print(f"Metric type: {metric_type.title()}")
    
    available_scoring_types = df['scoring_type'].dropna().unique()
    
    for scoring_type in ['pairs', 'dots', 'dots_only']:
        if scoring_type not in available_scoring_types:
            continue
            
        print(f"\n{scoring_type.upper().replace('_', ' ')} SCORING:")
        print("-" * 40)
        
        scoring_data = df[df['scoring_type'] == scoring_type]
        
        for base_label in base_labels:
            label_data = scoring_data[scoring_data['base_label'] == base_label]
            
            if len(label_data) == 0:
                continue
                
            print(f"\n{base_label}:")
            print(f"  Number of experiments: {len(label_data)}")
            
            # Find best performing temperature level
            best_temp = None
            best_mean = 0
            for temp in temp_levels:
                values = label_data[str(temp)].dropna()
                if len(values) > 0:
                    mean_val = values.mean()
                    if mean_val > best_mean:
                        best_mean = mean_val
                        best_temp = temp
            
            if best_temp:
                # Format output based on metric type or data inspection
                if metric_type and metric_type.lower() == 'percentage':
                    print(f"  Best temperature level: {best_temp} ({metric_type.lower()}: {best_mean:.1f}%)")
                elif metric_type:
                    print(f"  Best temperature level: {best_temp} ({metric_type.lower()}: {best_mean:.4f})")
                else:
                    # Auto-detect format
                    sample_values = []
                    for temp in temp_levels[:3]:
                        values = label_data[str(temp)].dropna()
                        sample_values.extend(values.tolist())
                    
                    if sample_values and max(sample_values) > 1:
                        print(f"  Best temperature level: {best_temp} (mean: {best_mean:.1f}%)")
                    else:
                        print(f"  Best temperature level: {best_temp} (mean: {best_mean:.4f})")
            
            # Calculate overall statistics across all temperature levels
            all_values = []
            for temp in temp_levels:
                values = label_data[str(temp)].dropna()
                all_values.extend(values.tolist())
            
            if all_values:
                if metric_type and metric_type.lower() == 'percentage':
                    print(f"  Overall mean: {np.mean(all_values):.1f}%")
                    print(f"  Overall std: {np.std(all_values):.1f}%")
                elif metric_type:
                    print(f"  Overall mean: {np.mean(all_values):.4f}")
                    print(f"  Overall std: {np.std(all_values):.4f}")
                else:
                    # Auto-detect format
                    if max(all_values) > 1:
                        print(f"  Overall mean: {np.mean(all_values):.1f}%")
                        print(f"  Overall std: {np.std(all_values):.1f}%")
                    else:
                        print(f"  Overall mean: {np.mean(all_values):.4f}")
                        print(f"  Overall std: {np.std(all_values):.4f}")

def main():
    parser = argparse.ArgumentParser(
        description='Plot temperature analysis with three scoring methods',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (auto-detect metric from data)
  python plot_scores_comparison.py results.csv
  
  # Specify metric type for y-axis labeling
  python plot_scores_comparison.py results_mean.csv --metric mean
  python plot_scores_comparison.py results_median.csv --metric median
  python plot_scores_comparison.py results_percentage.csv --metric percentage
  
  # Custom metric name
  python plot_scores_comparison.py results.csv --metric "Custom Score"
        """
    )
    
    parser.add_argument('csv_file', 
                        help='CSV file containing the results data')
    parser.add_argument('--metric', '-m', 
                        help='Metric type for y-axis labeling (e.g., "mean", "median", "percentage", "accuracy")')
    
    args = parser.parse_args()
    
    # Check if CSV file exists
    import os
    if not os.path.exists(args.csv_file):
        print(f"Error: CSV file '{args.csv_file}' not found")
        return 1
    
    # Run the plotting function
    try:
        plot_temperature_analysis(args.csv_file, args.metric)
        return 0
    except Exception as e:
        print(f"Error plotting data: {e}")
        return 1

if __name__ == "__main__":
    import sys
    
    # Check if running from command line with arguments
    if len(sys.argv) > 1:
        sys.exit(main())
    else:
        # Fallback to default behavior for backwards compatibility
        csv_file = "2n3r.csv"
        plot_temperature_analysis(csv_file)