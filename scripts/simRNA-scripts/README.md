# simRNA Tools

A comprehensive toolkit for [SimRNA](https://genesilico.pl/software/stand-alone/simrna/) REMC simulation management, structure analysis, and post-processing workflows for RNA 3D structure prediction.

---

## Installation

```bash
git clone <this-repo-url>
cd simRNA
```

Requirements:
- Python 3.6+
- SimRNA executable and data files
- Required Python packages: `numpy`, `pandas`, `matplotlib`, `seaborn`

---

## Overview

This repository provides tools for:

1. **Temperature Sweep Setup**: Automated preparation of SimRNA temperature sweep simulations
2. **Trajectory Processing**: Frame extraction and PDB structure generation
3. **Structure Analysis**: Secondary structure comparison and accuracy scoring
4. **Cluster Analysis**: Major cluster extraction and representative structure identification
5. **Batch Processing**: Parallel execution and result summarization

---

## Core Scripts

### 1. Temperature Sweep Setup

**`Temp_sweep.py`** - Creates folder structures for SimRNA temperature sweep simulations.

```bash
python Temp_sweep.py --output-dir temp_sweep_test \
                     --sequence-file sequences/3WJ_0_1_3 \
                     --structure-file sequences/3WJ_0_1_3.str \
                     --config-template templates/sim_config_template.in \
                     --temp-start 1.20 --temp-end 2.00 --temp-step 0.05 \
                     --create-run-script
```

**What it does:**
- Creates temperature-specific folders (e.g., `temp_1.20/`, `temp_1.25/`)
- Copies sequence and structure files to each folder
- Generates temperature-specific configuration files
- Creates symlinks to SimRNA data directory
- Optionally generates batch execution scripts

**Options:**
- `--temp-start`, `--temp-end`, `--temp-step`: Temperature range and step size
- `--simrna-data-path`: Path to SimRNA data directory for symlinks
- `--create-run-script`: Generate automated execution script

### 2. Frame Extraction and Processing

**`extract_low_temp_frames.py`** - Extracts frames from specific temperature ranges.

```bash
python extract_low_temp_frames.py simrna.log --base-name structure --min-temp 1 --max-temp 5 --output low_temp.trafl
```

**Options:**
- `--min-temp`, `--max-temp`: Temperature level range to extract (default: 0–10)
- `--min-frame`: Minimum frame number (default: 0)
- `--output`: Output file for extracted frames (default: `low_temp.trafl`)
- `--processes`: Number of parallel processes (default: auto)
- `--single-threaded`: Disable multiprocessing

**`extract_all_pdbs.py`** - Converts trajectory files to PDB structures.

```bash
python extract_all_pdbs.py /path/to/simrna/run --name structure_name
```

**What it does:**
- Runs `SimRNA_trafl2pdbs` for all trajectory files in directory
- Creates organized `pdbs/` directory
- Moves all generated PDB and `.ss_detected` files
- Provides detailed processing reports

**Options:**
- `--dry-run`: Show commands without executing
- `--suffixes`: Specific trajectory suffixes to process (e.g., "01,02,03")

### SASA calculation (sasa.py)

`sasa.py` computes solvent-accessible surface area (SASA) for selected nucleotide(s) using SimRNA trajectories and a reference PDB.

Basic usage:

```bash
python scripts/sasa.py --pdb reference.pdb --trafl traj.trafl --nucleotide C:1 --aa --probe-radius 1.4 --out sasa.csv --plot sasa.png --frames 1:500
```

Key options:
- `--pdb` : reference PDB used for mapping and AA reconstruction.
- `--trafl` : SimRNA `.trafl` trajectory file.
- `--nucleotide` : nucleotide(s) to analyze. Format: `A:10` (chain A resi 10) or `10` (1-based index). Multiple values comma-separated: `A:10,B:12`.
- `--aa` : attempt all-atom reconstruction via `SimRNA_trafl2pdbs` (requires SimRNA and `data/` available).
- `--frames` : frame specification passed to `SimRNA_trafl2pdbs` (e.g. `1:500` or `:` for all frames).
- `--out` : per-frame CSV output path (default `sasa_per_frame.csv`).
- `--plot` : write boxplots (absolute + relative SASA) to PNG.
- `--keep-pdb` : retain PDB files produced by the conversion (default: remove them after processing).

Notes:
- When `--aa` is used, `sasa.py` converts the selected frames to all-atom PDBs and computes SASA with FreeSASA; it will only process `_AA.pdb` files and will clean up intermediate coarse-grained PDBs and `.ss_detected` files.
- The script now calculates relative SASA (relative to an isolated residue reference) and writes a `rel_sasa` column to the CSV; plots include both absolute (Å^2) and relative SASA panels.
- If FreeSASA Python bindings are not installed, `--aa` will fail; install via `pip install freesasa`.

### SASA Summary Plotting (plot_sasa_summary.py)

`plot_sasa_summary.py` creates comparative visualizations from SASA boxplot summary data collected across multiple directories/runs.

Basic usage:

```bash
python scripts/plot_sasa_summary.py sasa_summary.csv --output comparison.png
```

Plot types and options:
- `--plot-type boxplot` : Box plots comparing median/quartiles across runs or nucleotides
- `--plot-type violin` : Violin plots showing distribution shapes (synthesized from summary stats)
- `--plot-type heatmap` : Heatmap of a specific metric across runs and nucleotides
- `--group-by run_id` : Group by run ID (compare nucleotides within each run)
- `--group-by nucleotide` : Group by nucleotide (compare runs for each nucleotide)
- `--metric abs|rel` : Plot absolute or relative SASA metrics
- `--nucleotides C:1 A:10` : Filter to specific nucleotides only
- `--create-all` : Generate all plot types plus summary table

Examples:

```bash
# Compare relative SASA across runs for each nucleotide
python scripts/plot_sasa_summary.py results.csv --metric rel --group-by nucleotide --plot-type boxplot

# Create heatmap of median absolute SASA
python scripts/plot_sasa_summary.py results.csv --plot-type heatmap --heatmap-metric median_abs

# Generate all plots and summary table
python scripts/plot_sasa_summary.py results.csv --create-all --output sasa_analysis
```

The script expects a CSV file with columns like `nucleotide`, `run_id`, `median_abs`, `mean_abs`, `std_abs`, `median_rel`, etc., as produced by `sasa.py --boxplot-csv`.

### 3. Structure Comparison and Scoring

**`compare_ss.py`** - Evaluates secondary structure prediction accuracy.

```bash
python compare_ss.py reference.str ss_detected_folder/ --output comparison.png --csv results.csv
```

**Scoring Modes:**
- **Base pairs only** (default): Score only Watson-Crick base pairs
- **With dots** (`--dots`): Include unpaired positions with half weight
- **Dots only** (`--dots-only`): Score only unpaired positions

**What it does:**
- Compares dot-bracket notation between reference and predicted structures
- Generates visualization plots showing accuracy across trajectories
- Exports detailed results with file names, scores, and accuracy ratios
- Groups files by trajectory for organized analysis

**Output Options:**
- `--all`: Create combined plot ordered by accuracy
- `--csv`: Export results to CSV format
- Multiple plot formats: individual trajectories or combined visualization

### 4. Temperature Analysis and Scoring

**`Analyze_scorings.py`** - Comprehensive temperature level analysis with all three scoring methods.

```bash
python Analyze_scorings.py --name structure --min-temp 1 --max-temp 16 --structure reference.str --summary results.csv --label experiment_1
```

**What it does:**
- Processes multiple temperature levels automatically
- Extracts frames, generates PDBs, and compares structures for each level
- Runs all three scoring methods (pairs, dots, dots_only) simultaneously
- Creates violin plots and comparison visualizations
- Updates summary CSV with results from all scoring methods

**Key Features:**
- **Automated workflow**: Complete pipeline from log files to final analysis
- **Three scoring methods**: Comprehensive structure evaluation
- **Continue on error**: Process remaining levels if one fails
- **Summary integration**: Automated CSV updates for batch analysis
- **Rich visualizations**: Violin plots and PDF reports

**Options:**
- `--continue-on-error`: Continue processing if individual levels fail
- `--skip-plots`: Skip plot generation for faster processing
- `--threshold`: Accuracy threshold for summary statistics (default: 0.8)
- `--median`, `--mean`: Use median/mean values instead of percentage above threshold
- `--quiet`: Reduce output verbosity

Position-specific analysis
--------------------------
`Analyze_scorings.py` and `compare_ss.py` support position-specific analysis. Use `--positions` (or `-p`) to restrict scoring/plots to specific nucleotide indices. Supported formats include:

- Single positions or comma-separated: `--positions 10,15,20`
- Ranges: `--positions 24-29`
- Mixed: `--positions 24-29,35,40-42`

When positions are provided, the tools will create position-specific trafl/pdb outputs (filename suffix `_pos_<ranges>`) and compute position-only accuracy statistics and plots. Example:

```bash
python scripts/Analyze_scorings.py --name structure --min-temp 1 --max-temp 16 --structure reference.str --positions "24-29" --summary results.csv --label experiment_positions
```

The `--positions` option is useful when you want accuracy metrics only for a sub-region of the molecule (e.g., a binding site or junction).

### 5. Results Summary and Visualization

**`create_summary.py`** - Creates summary files from existing accuracy CSV files.

```bash
python create_summary.py --summary results.csv --label experiment_1 --levels 1 2 3 4 5
```

**Features:**
- **File locking**: Thread-safe CSV updates for parallel processing
- **Auto-detection**: Automatically detect temperature levels from existing files
- **Multiple metrics**: Percentage above threshold, median, or mean accuracy
- **Three scoring methods**: Processes pairs, dots, and dots_only results

**Options:**
- `--threshold`: Accuracy threshold (default: 0.8)
- `--median`, `--mean`: Alternative summary metrics
- `--directory`: Custom directory containing CSV files

**`plot_scores_comparison.py`** - Advanced visualization of temperature analysis results.

```bash
python plot_scores_comparison.py results.csv --metric accuracy
```

**What it does:**
- Creates three-panel plots for pairs/dots/dots_only comparison
- Shows temperature trends with min-max ranges
- Generates publication-ready figures with legends
- Supports different metric types (accuracy, percentage, etc.)

### 6. Batch Processing and Automation

**`run_multiple_create_summary.sh`** - Batch processing for multiple result directories.

```bash
bash run_multiple_create_summary.sh _suffix results.csv --median
```

**Process Trajectories:**

```bash
bash process_trajectories.sh
```

**Extract Frames from Multiple Directories:**

```bash
bash extract_frames.sh
```

---

## Complete Workflow Examples

### 1. Temperature Sweep Analysis

```bash
# 1. Create temperature sweep folders
python Temp_sweep.py --output-dir temp_sweep \
                     --sequence-file RNA_sequence \
                     --structure-file RNA_sequence.str \
                     --config-template sim_config.in \
                     --temp-start 1.20 --temp-end 2.00 --temp-step 0.05 \
                     --create-run-script

# 2. Run SimRNA simulations (in each temperature folder)
cd temp_sweep && ./run_all_temps.sh

# 3. Analyze results across temperature range
python Analyze_scorings.py --name RNA_sequence \
                          --min-temp 1 --max-temp 16 \
                          --structure RNA_sequence.str \
                          --summary temperature_results.csv \
                          --label RNA_design_1
```

### 2. Structure Analysis Pipeline

```bash
# 1. Extract all PDB structures from trajectories
python extract_all_pdbs.py /path/to/simrna/run --name structure_name

# 2. Compare secondary structures with multiple scoring methods
python compare_ss.py reference.str pdbs/ --output comparison_pairs.png --csv results_pairs.csv
python compare_ss.py reference.str pdbs/ --output comparison_dots.png --csv results_dots.csv --dots
python compare_ss.py reference.str pdbs/ --output comparison_dots_only.png --csv results_dots_only.csv --dots-only

# 3. Create summary from results
python create_summary.py --summary final_results.csv --label structure_analysis --levels 1 2 3 4 5

# 4. Generate comparison plots
python plot_scores_comparison.py final_results.csv --metric accuracy

# 5. Calculate distances between nucleotides

python distance.py         --pdb level3_temp_frames-000001.pdb --trafl level3_temp_frames.trafl         --nt1 B:27 --nt2 C:2                  --out per_frame.csv --plot distances_2.png --boxplot-csv ../distance_2.csv --run-id Tr --aa > distance.log 2>&1 &
```

### 3. Cluster Analysis Workflow

```bash
# 1. Extract low temperature frames
python extract_low_temp_frames.py simulation.log --base-name structure --min-temp 1 --max-temp 5 --output low_temp.trafl

# 2. Perform clustering (using SimRNA tools)
clustering low_temp.trafl 0.01 15 >& clustering.log &

# 3. Extract major clusters
python extract_major_clusters.py -prefix low_temp_thrs15.00A -cutoff 5

# 4. Compare cluster representatives
python compare_ss.py reference.str cluster_structures/ --output cluster_comparison.png --csv cluster_results.csv --all
```

---

## Integration with RNA Pipeline

These simRNA tools integrate seamlessly with the main RNA Pipeline:

1. **Folder Preparation**: The RNA Pipeline's `simrna_converter.py` utility prepares simRNA-ready folders
2. **Structure Analysis**: Results from simRNA simulations can be analyzed using these tools
3. **Batch Processing**: Multiple pipeline results can be processed in parallel

**Example Integration:**

```bash
# From RNA Pipeline - prepare simRNA folders
python pipeline/utils/core/simrna_converter.py --config simrna_batch.json

# Transfer to compute cluster and run simulations
scp -i .pem -r simRNA/runs/* user@cluster:/path/to/simrna/

# After simulations complete - analyze results
python ../../../scripts/Analyze_scorings.py --name C2 --min-temp 1 --max-temp 16 --structure C2.str   --summary ../weights_tests.csv --label C2 --mean
