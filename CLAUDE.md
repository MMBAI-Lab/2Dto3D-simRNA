# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

This is **not** a SimRNA source tree. It is a workspace for running SimRNA (Genesilico, v3.20) simulations to go from RNA sequence / 2D structure to 3D ensembles, then analyzing breathing / conformational diversity of the resulting trajectories.

**Tracked in git:**

- [README.md](README.md) / [README_ES.md](README_ES.md) — parallel EN/ES user docs for the pipeline.
- [inputs/example/](inputs/example/) and [results/example/](results/example/) — a single reference run (inputs + post-processed outputs, no `.trafl`). Everything else under `inputs/` and `results/` is gitignored.

**Present on the local filesystem but gitignored:**

- `SimRNA_64bitIntel_Linux/` — upstream Linux binary distribution (ELF binaries + `data/` statistical potentials + user manual + academic license PDF). Cannot be redistributed publicly, so it is **not** in the repo; every user must download it from the official GeneSilico SimRNA site. Locally, treat its files as read-only.
- `docs/` — internal method report (`PRO - Method report_ simRNA.docx`) and a copy of the SimRNA user manual. Lab-internal, do not commit. The docx is the source of truth for protocol choices (16-replica REMC, 0.9–1.65 temperature range for TE-like structures, 15 Å clustering cutoff, `ETA_THETA_WEIGHT=0.4`). Read it before proposing parameter changes.

## Binaries and how to invoke them

All binaries require a `data/` directory (or symlink) in the **current working directory** — they will refuse to start otherwise. The canonical setup is a symlink:

```
ln -s /path/to/SimRNA_64bitIntel_Linux/data data
```

Binaries (in [SimRNA_64bitIntel_Linux/](SimRNA_64bitIntel_Linux/)):

- `SimRNA` — main Monte Carlo engine. Inputs: `-s <seq>` **or** `-p <pdb>` (mutually exclusive; with `-p`, sequence comes from the PDB). Optional: `-S <ss_restraints>`, `-c <config>`, `-r <restraints>`, `-o <prefix>`, `-n <iterations>`, `-R <seed>`, `-E <nReplicas>` to enable Replica Exchange, `-P` to treat PDB occupancy/bfactor as per-atom position restraints, `-l` to additionally emit per-frame PDBs.
- `clustering file.trafl <fraction> <rmsd_thrs1> [<rmsd_thrs2> ...]` — RMSD clustering on a trajectory. Produces `*_thrs<X>A_clust<N>.trafl` files.
- `SimRNA_trafl2pdbs <template.pdb> <file.trafl> [AA] <range>...` — extract frames to PDBs. Ranges: `22` (single), `103:153` (inclusive), `:22` (start→22), `1022:` (1022→end). The `AA` keyword reconstructs full-atom from the 5-pseudoatom coarse-grained representation.
- `calc_rmsd_to_1st_frame <file.trafl> <out.rmsd_e>` — RMSD of each frame to frame 1.
- `traflView` — trajectory viewer.
- `trafl_extract_lowestE_frame.py <file.trafl>` — **Python 2 only** (uses `print >>sys.stderr` syntax). Writes `<name>_minE.trafl`. Do not "port" it unless asked; the file is part of the vendored distribution.

## Typical run (REMC from sequence + dot-bracket restraints)

```
ln -s /path/to/SimRNA_install/data data
nohup SimRNA -s example -S example.str -c sim_config.in -E 16 -o example >& example.log &
```

A minimal `sim_config.in` ships with the SimRNA distribution (`config.dat` / `config_n0.dat`); see the docx for the parameter dictionary. Key outputs: `*.trafl` (trajectory), `*.pdb` (starting frame), `*.ss_detected` (detected 2D), and `example.log` (energies + replica exchange events).

## Downstream analysis scripts are NOT in this repo

The method report references a suite of Python scripts — `extract_low_temp_frames.py`, `extract_major_clusters.py`, `Analyze_scorings.py`, `plot_scores_comparison.py`, `extract_all_pdbs.py`, `compare_ss.py`, `multi-pdb.py` — that parse REMC logs, score 2D agreement, and build multi-model PDBs for ChimeraX. **None of them are checked in here.** If a task requires one, ask the user where the `simRNA/` helper tree lives before assuming paths. The docx shows invocation patterns but not the implementations.

## Conventions worth knowing before editing

- Temperatures are on a **relative** scale (T≈1 is standard). "Around T=1" in the protocol means REMC levels 2–4 when using the documented 16-level 0.9–1.65 ladder — the indices, not the temperatures, are what the analysis scripts key on.
- Broken `data` symlinks are the single most common failure mode for these binaries. If a simulation dies at startup, check the symlink before anything else.
- SS restraint weight (`SECOND_STRC_RESTRAINTS_WEIGHT`) is the main knob being tuned; the docx describes the too-low / too-high diagnostic signatures. Don't change it without reading that section.
- **Bilingual READMEs must stay in sync.** [README.md](README.md) (English) and [README_ES.md](README_ES.md) (Spanish) are parallel translations of the same content. Any edit to one must be mirrored in the other in the same change — never update only one language. If the user requests a change in one language, produce both edits in parallel.
- **Do not git-add `SimRNA_64bitIntel_Linux/` or `docs/`.** Both are in `.gitignore` and must stay that way — the SimRNA distribution has a non-redistribution license, and `docs/` is lab-internal. If `git status` ever shows those paths as untracked or modified, something broke in `.gitignore`.
- **Reference run is at `inputs/example/` + `results/example/` only.** Any other per-run directory under those is ignored by design; do not add whitelist rules for new ones. `*.trafl` is excluded everywhere — the repo stores post-processed outputs only.
