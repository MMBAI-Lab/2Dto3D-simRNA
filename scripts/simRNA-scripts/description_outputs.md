# Output Folder Organization and Naming Conventions

This document describes the organization and naming conventions used in the output folder.

---

## Naming Conventions

- **_ss{number}**: Indicates the second structure restraints weight used in the run.
  - `1` = 1.0
  - `05` = 0.5
  - `005` = 0.05
  - etc.
- **_rep{number}**: Identifies individual repetitions of a run (e.g., `_rep1`, `_rep2`, ...).

---

## Example

To find results for a 2n3r run with a second structure restraint weight of 0.05, repetition 2:

- Go to: `output/2n3r/2n3r_ss005_rep2/`

## Structure

The output folder is organized in different folders that include structures dirs like TE, BL and TR, controls like 2n3r, temp_sweeps used for parameter optimization and an others file with other runs.
