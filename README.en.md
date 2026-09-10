# Auto-MLWFs Skill

[简体中文](README.md) | **English**

[![Validate Skill](https://github.com/sqwang2/auto-mlwfs/actions/workflows/ci.yml/badge.svg)](https://github.com/sqwang2/auto-mlwfs/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/Python-3.11%2B-3776AB?logo=python&logoColor=white)](skills/auto-mlwfs/requirements.txt)
[![License: MIT](https://img.shields.io/badge/License-MIT-2A7F62)](LICENSE)

## Installation

With Node.js/npm and Git installed, run this command to install the skill for **Codex**:

```bash
npx skills add sqwang2/auto-mlwfs --skill auto-mlwfs -a codex -g -y
```

For **Claude Code**:

```bash
npx skills add sqwang2/auto-mlwfs --skill auto-mlwfs -a claude-code -g -y
```

`-g` installs the skill at the user level; omit it to install into the current project. After installation, ask your agent to use `auto-mlwfs`. Start a new session if the skill does not appear immediately.

<details>
<summary><strong>Stuck at Cloning repository? Install from GitHub's official source archive</strong></summary>

If your terminal cannot connect to `github.com` but can reach `codeload.github.com`, press `Ctrl+C` to stop the stalled command, then run:

```bash
npx --yes skills@latest add https://codeload.github.com/sqwang2/auto-mlwfs/tar.gz/refs/heads/main --skill auto-mlwfs -a codex -g -y
```

This downloads and installs the source archive directly from GitHub, bypassing Git cloning. It is not a third-party mirror; changing the npm registry will not resolve connectivity problems during the GitHub clone step. For Claude Code, replace `-a codex` with `-a claude-code`.

</details>

Installation uses the third-party [skills CLI](https://github.com/vercel-labs/skills) and includes the skill instructions, scripts, and configuration. Prepare Python dependencies and external programs such as VASP and Wannier90 separately, as described in [Requirements](#requirements).

## Overview

A VASP-Wannier90 workflow for maximally localized Wannier functions (MLWFs), designed for AI agents.

This project brings projection selection, VASP/Wannier90 input preparation, energy-window construction, fitting assessment, and Bayesian optimization into a reusable skill. It starts from a completed VASP SCF calculation and guides the subsequent Wannierization workflow, with the agent and user confirming compute resources and job submission authorization together.

> Current status: research and engineering validation. Test the workflow on an existing case before applying it to a new material system.

## Features

- Automatically recommend Wannier projections from DOS/PDOS in `vasprun.xml`.
- Prepare inputs for scalar, spin-orbit coupling (SOC), and noncollinear calculations.
- Prepare VASP band, projection, and Wannier90 inputs automatically.
- Set `NBANDS` from `NELECT`, `NUM_WANN`, and the actual MPI rank count.
- Construct initial disentanglement and frozen windows from the projection band energies.
- Generate DFT-Wannier fitting plots and machine-readable fitting assessments with LESA.
- Apply Bayesian optimization (BO) to the four energy-window boundaries when the initial result fails assessment.
- Let the agent adapt to the current Slurm or PBS environment instead of relying on fixed submission templates.
- Return structured results to the coordinating agent through `stage_summary.json` and BO batch JSON files.

## Workflow

```text
completed VASP SCF
        |
        +--> DOS/PDOS analysis --> projection recommendation
        |
        +--> VASP line-mode k-path preparation
        |
        +--> VASP projection calculation
        |
        +--> Wannier90 input and physics-window preparation
        |
        +--> Wannier90 calculation
        |
        +--> final VASP band calculation on Wannier90 k-points
        |
        +--> fitting plot + LESA assessment
                         |
                         +--> pass: deliver results
                         |
                         `--> fail: four-dimensional BO
```

## Repository Structure

```text
auto-mlwfs/
├── README.md                        # Chinese documentation
├── README.en.md                     # English documentation
├── LICENSE
├── CITATION.cff
├── skills/
│   └── auto-mlwfs/
│       ├── SKILL.md                 # Core agent instructions
│       ├── LICENSE
│       ├── requirements.txt
│       ├── configuration/
│       │   ├── projection.yaml      # Projection selection parameters
│       │   ├── fitting.yaml         # LESA fitting criteria
│       │   └── optimization.yaml    # Default BO configuration
│       ├── docs/                    # Workflows, criteria, and troubleshooting
│       ├── examples/agi/
│       │   └── optimization_phase2.yaml
│       └── scripts/
│           ├── dos_and_analysis.py
│           ├── vasp_band.py
│           ├── prepare_projection.py
│           ├── prepare_wannier.py
│           ├── prepare_band_from_wannier_kpt.py
│           ├── plot_fitting.py
│           ├── lesa_fitting.py
│           ├── optimize_windows.py
│           ├── check_success.py
│           └── ...                 # Shared parsing and output helpers
└── tests/                           # Installation portability and synthetic-data tests
```

## Requirements

### External Programs

- VASP with an appropriate license.
- VASPKIT for automatic high-symmetry k-path generation; it can be skipped when a valid path is already available.
- Wannier90.
- An optional Slurm or PBS scheduler environment.

Users supply program locations, launch commands, queues, core counts, memory, and environment modules for their execution environment. The project does not embed site-specific paths.

VASP binaries, POTCAR files, other external programs, and cluster credentials are not distributed with this project.

### Python Dependencies

Use **Python 3.11 or newer**. Basic analysis and plotting require:

```text
numpy >= 2.0
pymatgen
matplotlib
PyYAML
```

Bayesian optimization additionally requires:

```text
scikit-optimize
```

For a user-level Codex installation through the skills CLI on macOS/Linux, create an isolated Python environment:

```bash
export AUTO_MLWFS_SKILL_DIR="$HOME/.agents/skills/auto-mlwfs"
python3 -m venv "$HOME/.venvs/auto-mlwfs"
source "$HOME/.venvs/auto-mlwfs/bin/activate"
python -m pip install -r "$AUTO_MLWFS_SKILL_DIR/requirements.txt"
```

Confirm that `python3` is Python 3.11+. Other agents, project-level installations, or older installers may use different directories. Set `AUTO_MLWFS_SKILL_DIR` to the directory reported by the installer that contains `SKILL.md`.

See [`requirements.txt`](skills/auto-mlwfs/requirements.txt) for the full dependency version ranges. When running on a cluster, both the skill directory and Python environment must be accessible from the execution nodes; job scripts should use the corresponding interpreter.

## Input Requirements

A typical project starts from a completed `scf/` directory. Depending on the stage, the following files are usually required:

```text
scf/
├── INCAR
├── POSCAR
├── KPOINTS
├── POTCAR                           # Supplied by the user under its license
├── OUTCAR
├── vasprun.xml
├── EIGENVAL
└── CHGCAR
```

DOS/PDOS analysis requires `vasprun.xml` to contain projected density-of-states data.

The skill installation directory is independent of the calculation directory. Run the commands below from the calculation project root, using `AUTO_MLWFS_SKILL_DIR` defined during environment setup:

```text
project/
└── scf/
```

For example, ask your agent: "Use auto-mlwfs to analyze the PDOS in this project's scf/ directory and recommend projection orbitals. Check the inputs and environment before preparing subsequent calculations."

## Quick Start

### 1. Analyze DOS/PDOS and Recommend Projections

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/dos_and_analysis.py" \
  --workdir scf \
  --config "$AUTO_MLWFS_SKILL_DIR/configuration/projection.yaml"
```

Main outputs in `scf/outputs/`:

- `projection_summary.json`
- `projections.txt`
- `dos_projected.csv`
- `dos_projected_weights.csv`
- `DOS.png`

### 2. Prepare the High-Symmetry k Path

First, ask the user for the actual MPI rank count that will run VASP:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/vasp_band.py" --root . --run-np <actual-mpi-ranks>
```

This stage uses VASPKIT to generate a line-mode k path for Wannier90. It is not the final VASP band calculation used for fitting.

### 3. Prepare and Run the Projection Calculation

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_projection.py" \
  --root . \
  --run-np <actual-mpi-ranks>
```

The script reads `scf/outputs/projections.txt`, computes `NUM_WANN`, and prepares `proj/`. The agent then detects Slurm/PBS, asks for missing execution parameters, creates the submission script, and submits only after explicit user confirmation.

Check the completed calculation:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/check_success.py" --stage projection --workdir proj
```

### 4. Prepare and Run Wannier90

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_wannier.py" --root .
```

The script constructs the initial energy windows from the projection band energies, writes `wann/wannier90.win`, and generates `wann/window_summary.json`. Run Wannier90 using the current cluster environment.

Check the completed calculation:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/check_success.py" --stage wannier --workdir wann
```

### 5. Prepare the Final VASP Band Calculation

After Wannier90 produces `wannier90_band.kpt`, run:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_band_from_wannier_kpt.py" \
  --root . \
  --wann-dir wann \
  --band-dir band \
  --run-np <actual-mpi-ranks>
```

The script rewrites `band/KPOINTS` using the k points actually sampled by Wannier90. The agent and user then confirm resources and submission authorization before running the final VASP band job.

Check the completed calculation:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/check_success.py" --stage vasp-band --workdir band
```

### 6. Plot the Bands and Run LESA Assessment

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/plot_fitting.py" --root . --vasp-dir band --wann-dir wann

python "$AUTO_MLWFS_SKILL_DIR/scripts/lesa_fitting.py" \
  --root . \
  --config "$AUTO_MLWFS_SKILL_DIR/configuration/fitting.yaml"
```

Main outputs:

- `<project-name>_fitting.png`
- `outputs/fitting_lesa_summary.json`
- `outputs/fitting_lesa_report.txt`

The fitting plot supports human inspection; the LESA JSON provides the agent's default machine-readable assessment. Retain both together.

The current fitting parser cannot fully assess both spin channels in collinear spin-polarized `ISPIN=2` calculations. Input preparation support therefore does not imply end-to-end validation for every spin mode.

## NBANDS Rules

`NBANDS` is recalculated from `NELECT` in `scf/OUTCAR` rather than inherited directly from the previous INCAR.

For the band and final-band stages:

```text
NBANDS = ceil_to_multiple(2 × NELECT, actual_mpi_ranks)
```

For the projection stage:

```text
NBANDS = ceil_to_multiple(max(2 × NELECT, NUM_WANN), actual_mpi_ranks)
```

These rules use the MPI ranks that actually run the program, which may differ from the total core count requested from the scheduler. If memory constraints require changing the actual rank count, rerun the corresponding preparation script.

## Scheduler and Submission Authorization

The project does not require a cluster-specific submission YAML. The agent should:

1. Detect whether the environment uses Slurm or PBS and consult existing successful job scripts.
2. Ask the user for VASP/Wannier90 commands, queue, requested cores, actual MPI ranks, memory, wall time, and environment setup.
3. Determine from INCAR whether to use standard VASP or the SOC/noncollinear executable.
4. Present the resolved resources and launch command.
5. Submit only after explicit user confirmation.

If a band or projection job shows credible evidence of insufficient memory, the agent asks the user whether to reduce actual MPI ranks, request more memory, or do both. It does not resubmit automatically.

## Directory Collision Protection

For newly created `band/`, `proj/`, and `wann/` directories:

- If the target directory does not exist, use the requested name.
- If the target directory exists, write the new results to `<name>_agent/`.
- If both the original directory and `<name>_agent/` exist, stop and ask the user.

Subsequent commands should read the actual `workdir` from `stage_summary.json` and pass it through `--band-dir`, `--proj-dir`, or `--wann-dir` to avoid overwriting existing results.

## Bayesian Optimization

When the initial energy windows fail LESA assessment, the following boundaries can be optimized:

- `dis_win_min`
- `dis_win_max`
- `dis_froz_min`
- `dis_froz_max`

Before generating each new batch, the agent must ask how many candidates the user wants to run concurrently. The recommendation is 10, but the actual count follows the user's choice:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" \
  --root . \
  propose --batch-size <user-choice>
```

The command creates scheduler-independent trial directories and `optimization/batch-XXX.json`. A scheduler job runs each candidate with:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" \
  --root . \
  run-trial \
  --trial-id <trial-id> \
  --wannier-command '<launcher> <wannier90.x> {{ seedname }}'
```

After all candidates finish:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" --root . collect
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" --root . status
```

`skills/auto-mlwfs/examples/agi/optimization_phase2.yaml` is a second-phase configuration created after `trial-0006` reached the search boundary during the first AgI phase. It is not a general default and should not be applied directly to other materials.

The default `bo.max_evaluations` is `null`, which imposes no evaluation limit. Set an explicit limit in your project configuration when you need to control the compute budget. When using a custom `--config`, use the same configuration for `propose`, `run-trial`, `collect`, and `status`.

## Machine-Readable Outputs

The main preparation scripts write `stage_summary.json` in the actual working directory, recording:

- The current stage and status.
- The actual working directory and whether the `_agent` suffix was used.
- Key parameters such as `NBANDS`, `NUM_WANN`, VASP mode, and initial energy windows.
- Scheduler inputs that still require user confirmation.
- The recommended next action.

The BO file `optimization/batch-XXX.json` also records candidate parameters, trial IDs, working directories, and constraint counts. Coordinating agents should read these JSON files in preference to parsing console text.

## Configuration and Detailed Rules

Bundled configurations are located relative to the installed skill directory. To customize parameters, copy the YAML into your calculation project and select it with `--config`:

- LESA and BO resolve relative `--config` paths against `--root`; DOS resolves relative `--config` paths against the current shell directory.
- Calculation data paths inside YAML are relative to `--root`, or to `--workdir` for DOS.
- BO uses the bundled LESA configuration when `paths.fitting_config: null`. An explicit configuration path is relative to `--root`, unless it is absolute.

Detailed guidance for each stage:

- Projection rules: [`docs/projection_rules.md`](skills/auto-mlwfs/docs/projection_rules.md)
- Band workflow: [`docs/band_workflow.md`](skills/auto-mlwfs/docs/band_workflow.md)
- Scheduler interaction: [`docs/scheduler_guidance.md`](skills/auto-mlwfs/docs/scheduler_guidance.md)
- LESA criteria: [`docs/lesa_fitting.md`](skills/auto-mlwfs/docs/lesa_fitting.md)
- Success criteria: [`docs/success_criteria.md`](skills/auto-mlwfs/docs/success_criteria.md)
- Troubleshooting: [`docs/troubleshooting.md`](skills/auto-mlwfs/docs/troubleshooting.md)

## Usage Notes and Feedback

This initial skill was developed from practical materials calculation workflows. Different clusters, VASP/Wannier90 builds, and material systems may reveal further compatibility issues. Preserve inputs, logs, `stage_summary.json`, LESA outputs, and the exact commands that trigger a problem to support reproduction and improvement.

Submit problems and suggestions through [GitHub Issues](https://github.com/sqwang2/auto-mlwfs/issues). Remove cluster credentials and license-restricted files before sharing a minimal reproduction.

## Development and Validation

With dependencies installed, run these commands from the repository root:

```bash
python -m unittest discover -s tests -v
npx skills add . --list
```

The existing tests cover skill installation portability and analysis behavior on synthetic band data. CI runs on Python 3.11 and 3.12. Passing tests does not establish scientific validation across all materials, external program versions, or clusters.

## License and Citation

This project uses the [MIT License](LICENSE). External scientific software retains its own licensing requirements.

When using this software in research, cite the project using [`CITATION.cff`](CITATION.cff), and cite the original literature for VASP, Wannier90, and the methods you use as appropriate.

---

[简体中文](README.md) | **English**
