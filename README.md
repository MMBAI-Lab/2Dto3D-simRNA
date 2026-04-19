# 2Dto3D-simRNA

End-to-end pipeline to go from an RNA **sequence + dot-bracket secondary structure** to an analyzed **3D conformational ensemble** using SimRNA with Temperature Replica Exchange Monte Carlo (TREMD / REMC) sampling.

A Spanish version of this document is available at [README_ES.md](README_ES.md).

## Goal

Given:
- an RNA sequence (one or more chains, standard A/C/G/U codes), and
- a predicted secondary structure in multi-line dot-bracket notation,

produce:
- a REMC trajectory at multiple temperatures around the standard relative `T ≈ 1`,
- a set of representative 3D conformations (cluster centroids),
- quantitative scoring of how well each sampled ensemble reproduces the target 2D structure,
- diagnostic plots to rank design variants by stability, diversity, and conservation of paired / unpaired regions.

The target use cases are (a) **3D validation** of RNA designs (do the proposed sequences actually fold as intended under the chosen restraint regime?) and (b) **breathing analysis** (how much does the 3D structure fluctuate around the designed state at physiological-like temperatures?).

## Why SimRNA

SimRNA is a coarse-grained Monte Carlo engine developed by GeneSilico that represents each nucleotide with five pseudoatoms (two backbone: P, C4'; three base orientation). It uses statistical potentials derived from solved RNA structures and supports single- and multi-chain RNAs. Solvent is implicit — not vacuum, but not explicit water either, which can introduce artifacts and motivates the analysis protocol in this repo.

REMC sampling runs several replicas in parallel at different relative temperatures and periodically exchanges them to escape local minima. The trajectory at the lowest temperatures captures the designed fold; higher-temperature replicas probe the unfolding landscape, which lets us check both that the design is stable near `T = 1` and that it actually unfolds when thermally stressed.

## Repository layout

```
.
├── inputs/
│   └── example/               # reference run inputs (tracked)
├── results/
│   └── example/               # reference run post-processed results (tracked, no .trafl)
├── CLAUDE.md                  # notes for Claude Code sessions on this repo
├── README.md                  # this file
└── README_ES.md               # Spanish version
```

Only the single reference run under `inputs/example/` and `results/example/` is tracked. Any other directory you create under `inputs/` or `results/` is gitignored by default — use them freely for your own runs without polluting the repo. The heavy trajectory files (`*.trafl`) are excluded even for the reference example; the repo keeps only the light-weight post-processed outputs (scoring CSVs, centroid PDBs, plots).

The **SimRNA academic distribution** (binaries + `data/` statistical potentials) is **not** distributed from this repo — the license does not permit public redistribution. Download it separately before running the pipeline; see [Installing SimRNA](#installing-simrna) below.

## Installing SimRNA

Download the SimRNA academic distribution (v3.20 for 64-bit Linux is what this pipeline is calibrated against) from the official GeneSilico / IIMCB SimRNA website. The distribution ships:

- `SimRNA` — main MC / REMC engine
- `clustering` — RMSD clustering of a trajectory
- `SimRNA_trafl2pdbs` — extract frames from a `.trafl` to PDB (optionally with full-atom reconstruction via `AA`)
- `calc_rmsd_to_1st_frame` — per-frame RMSD to the first frame
- `traflView` — trajectory viewer
- `trafl_extract_lowestE_frame.py` — **Python 2** helper that pulls out the lowest-energy frame
- `data/` — statistical-potential histograms required at runtime
- `SimRNA_academic_license.pdf` — read before use; redistribution is restricted

After download, either add the distribution directory to your `PATH` or invoke the binaries by absolute path. Wherever you run a simulation from, SimRNA needs a `data/` directory (or symlink) in the current working directory — see [Stage 1](#stage-1--prepare-the-run-directory).

## Inputs

1. **Sequence file** — plain text, nucleotides `A/C/G/U`; multiple chains are separated by spaces.
2. **Secondary-structure restraints** (`*.str`, optional but used throughout this pipeline) — multi-line dot-bracket string matching the sequence, one line per chain.
3. **Configuration file** (`*.in` / `*.dat`) — SimRNA parameters. The SimRNA distribution ships a minimal `config.dat` / `config_n0.dat` example next to the binaries. Typical TREMD values:

   ```
   NUMBER_OF_ITERATIONS 16000000
   TRA_WRITE_IN_EVERY_N_ITERATIONS 16000
   INIT_TEMP 1.65
   FINAL_TEMP 0.90
   BONDS_WEIGHT 1.0
   ANGLES_WEIGHT 1.0
   TORS_ANGLES_WEIGHT 0.0      # deprecated, keep at 0
   ETA_THETA_WEIGHT 0.4
   SECOND_STRC_RESTRAINTS_WEIGHT 0.3
   ```

4. **`data` symlink** — SimRNA refuses to start without a `data/` directory (or symlink) in the current working directory. This is the single most common source of run failures. Point it at the `data/` directory that shipped with your SimRNA download:

   ```bash
   ln -s /absolute/path/to/your/SimRNA_install/data data
   ```

5. **(Optional) starting PDB** — with `-p` the sequence is read from the PDB; with `-P` occupancy/B-factor are interpreted as per-atom position restraints.

## Parameter notes

- `NUMBER_OF_ITERATIONS` × `TRA_WRITE_IN_EVERY_N_ITERATIONS`: the ratio sets how many frames (and, in REMC, how many exchange attempts) you get. Common working points: 5M / 5k / 20 replicas or 16M / 16k / 10 replicas.
- `INIT_TEMP` / `FINAL_TEMP`: relative scale — `T ≈ 1` is the physiological operating point. For TE-like structures the protocol uses 16 replicas spanning roughly `0.9 → 1.65`. The high end has to actually unfold the structure or REMC is not doing its job.
- `SECOND_STRC_RESTRAINTS_WEIGHT`: main knob.
  - **Too low:** the ensemble near `T = 1` drifts away from the target 2D structure.
  - **Too high:** the structure survives at the highest temperatures (REMC stops working) and variability at `T ≈ 1` collapses (not useful for breathing analysis).

  Start around `0.05` for a single exploratory run, sweep upward if agreement is too weak, downward if the high-T replicas never unfold. Variability depends on this weight, so run **≥ 3 repetitions at the chosen weight** to estimate noise.
- `TORS_ANGLES_WEIGHT` is ignored in current SimRNA and has been superseded by `ETA_THETA_WEIGHT`; leave it at `0`.

## Pipeline

### Stage 1 — Prepare the run directory

```bash
mkdir -p results/my_run && cd results/my_run
ln -s /abs/path/to/your/SimRNA_install/data data
cp ../../inputs/my_run/example            .     # sequence
cp ../../inputs/my_run/example.str        .     # dot-bracket restraints
cp ../../inputs/my_run/sim_config.in      .
```

> Only `inputs/example/` and `results/example/` are tracked in git; everything else under those two directories is local-only (see [Repository layout](#repository-layout)).

### Stage 2 — Run REMC (TREMD)

Three invocations cover the typical cases; **(A) is the one used most often**.

**(A) REMC from sequence + 2D restraints — primary workflow (16 replicas):**

```bash
nohup SimRNA -s example -S example.str -c sim_config.in -E 16 -o example >& example.log &
```

**(B) MC from a starting PDB + 2D restraints (single-temperature):**

```bash
nohup SimRNA -p example.pdb -S example.str -c sim_config.in -o example >& example.log &
```

**(C) REMC from sequence, no restraints (20 replicas):**

```bash
nohup SimRNA -s example -c sim_config.in -E 20 -o example >& example.log &
```

Outputs:

- `example_<NN>.trafl` — trajectory per replica (coarse-grained, one frame per `TRA_WRITE_IN_EVERY_N_ITERATIONS` iterations)
- `example_<NN>-000001.pdb` — first-frame PDB per replica (used as template for `trafl2pdbs`)
- `example.ss_detected` — secondary structure detected in the initial conformation
- `example.log` — energy trace and replica-exchange record

### Stage 3 — Extract low-temperature frames around `T ≈ 1`

For the standard 16-level `0.9 → 1.65` ladder, the replicas that sit around `T = 1` are typically levels `2`, `3`, and `4` (relative temperatures `0.95`, `1.00`, `1.05`). Pool their frames into a single trajectory:

```bash
python simRNA/extract_low_temp_frames.py example.log \
  --base-name example --min-temp 2 --max-temp 4 --output around_1.trafl
```

### Stage 4 — Cluster

Cluster the pooled low-temperature frames with a 15 Å RMSD cutoff (adjust per system):

```bash
clustering around_1.trafl 1.0 15 >& clustering.log &
```

This produces `around_1_thrs15.00A_clust<NN>.trafl` files. Frames inside each cluster are ordered by similarity to the centroid; per-frame RMSDs to the centroid are in `clustering.log`.

### Stage 5 — Extract centroids / representative PDBs

```bash
# Centroids of all clusters above a representation cutoff
python extract_major_clusters.py --cutoff <percent>

# Or, the centroid (first frame) of a specific cluster, full-atom reconstructed:
SimRNA_trafl2pdbs example_01-000001.pdb around_1_thrs15.00A_clustXX.trafl 1 AA
```

### Stage 6 — Score 2D agreement

Compare the secondary structure observed in the ensemble against the target dot-bracket:

```bash
# Full per-temperature scoring + across-level comparison, appending a summary row
python /path/to/Analyze_scorings.py \
  --name example --min-temp 1 --max-temp 16 \
  --structure example.str --label example \
  --summary ../summary.csv --threshold 0.8
```

For each temperature this emits three scores: base-pair agreement, full-sequence agreement (paired + unpaired), and dots-only agreement (unpaired-only). The dots and all-sequence scores are **not expected to go to zero at high temperatures** — unpaired regions are what high-T samples are supposed to look like.

Compare multiple directories (e.g. replicates or design variants):

```bash
python /path/to/plot_scores_comparison.py summary.csv
```

Manual alternative for the same step:

```bash
python /path/to/extract_all_pdbs.py /path/to/dir --name example --suffixes 01,02,03,04,05
python /path/to/compare_ss.py example.str pdbs/
```

### Stage 7 — Inspect trajectories in ChimeraX or VMD

Combine per-frame PDBs into a multi-model PDB that both ChimeraX and VMD open as a trajectory:

```bash
python multi-pdb.py --folder ./pdbs --output trajectory.pdb
```

Open it in ChimeraX with `File > Open` (models step through automatically), or in VMD with `vmd trajectory.pdb` (each model becomes a frame).

## Diagnostic flow (what to look for)

1. **Single exploratory run** with a low restraint weight (`~0.05`), 16 replicas, `0.9 → 1.65`.
2. Run scoring. Near `T = 1` (levels 2–4) base-pair agreement should be high without being pinned to 1.0; at the two highest levels (15–16) the structure should be mostly lost.
3. If the target structure is not recovered near `T = 1` → **increase** `SECOND_STRC_RESTRAINTS_WEIGHT`.
4. If the structure survives at the highest temperatures or variability at `T = 1` collapses → **decrease** the weight.
5. Once a weight is chosen, repeat the run **≥ 3×** to quantify noise at that weight.
6. Run all design variants, do the rupture analysis and whisker plot over base-pairs / dots / base-pairs+dots. Pre-select devices on stability, diversity, and conservation of unpaired regions.
7. Cluster the selected variants (15 Å cutoff as starting point) and take centroids of clusters with representation above a chosen cutoff (e.g. 5%).
8. **(Optional)** re-run at a tuned temperature range for extra sensitivity.

## Caveats

- **Solvation is implicit.** The coarse-grained energy function does not see water explicitly; treat absolute energies as relative and validate conclusions with ensemble-level statistics rather than single frames.
- **Sequence length drives cost** super-linearly; keep chains as short as the biological question allows.
- **Broken `data` symlinks** are the top cause of silent restarts of this pipeline — check them first when a run refuses to start or crashes immediately.
- **The downstream Python analysis scripts** (`extract_low_temp_frames.py`, `extract_major_clusters.py`, `Analyze_scorings.py`, `plot_scores_comparison.py`, `extract_all_pdbs.py`, `compare_ss.py`, `multi-pdb.py`) are **not yet vendored in this repository** — they live in an external `simRNA/` helper tree and are referenced here by the documented protocol.

