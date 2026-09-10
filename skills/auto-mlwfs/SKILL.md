---
name: auto-mlwfs
description: Build automated VASP to Wannier90 MLWF workflows, especially selecting Wannier projection orbitals from SCF vasprun.xml DOS/PDOS results, generating projection recommendations, and preparing downstream wannier90 inputs.
---

# Auto MLWFs

Use this skill when a VASP SCF calculation must be analyzed before Wannierization.

## Installed Skill and Runtime

Locate this `SKILL.md` in the agent's installed skill directory and use its
containing directory as `AUTO_MLWFS_SKILL_DIR` in the shell examples below.
Set that variable to the actual absolute path before running commands; do not
assume the calculation contains a `skill/` directory. Keep calculation data
in the user's project directory, separate from the installed skill.

Use Python 3.11 or newer with the packages in `requirements.txt`. Install them
in an isolated environment when needed, then use that interpreter for all
stages and scheduler jobs. VASP, Wannier90, VASPKIT, MPI, and cluster access
are external prerequisites, configured for the user's environment.

The default YAML files are loaded from this skill's `configuration/` directory.
For LESA and optimization, explicit `--config` paths and calculation paths
inside YAML are relative to `--root` unless absolute. For DOS, explicit
`--config` paths are relative to the shell's current directory, and data paths
inside YAML are relative to `--workdir`. In optimization YAML,
`paths.fitting_config: null` selects the bundled fitting configuration.
Pass explicit overrides to reuse older project-local configuration files.

Copy configurations into the calculation project before customizing them;
do not modify the installed defaults. Preserve the chosen `--root`, `--config`,
runtime interpreter, and installed script path in subsequent scheduler jobs.

## DOS Projection Analysis

1. Work from a completed VASP SCF directory containing `vasprun.xml`.
2. Run `scripts/dos_and_analysis.py` with `configuration/projection.yaml`.
3. Use the generated files in the SCF `outputs/` directory:
   - `dos_projected.csv`: full element/orbital/spin-resolved PDOS table.
   - `dos_projected_weights.csv`: integrated band-edge weights by element and orbital.
   - `projection_summary.json`: machine-readable band-edge analysis and ranked channels.
   - `projections.txt`: one `Element:orbital` projection per line for `wannier90.win`.
   - `DOS.png`: selected useful PDOS channels plotted automatically.

Example:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/dos_and_analysis.py" --workdir .
```

## Projection Rule

Use `docs/projection_rules.md` for the current decision rule. In short:

- Find the connected nonzero-DOS manifold around `E_F`, expanding DOS segments
  by half the band gap plus the configured padding before merging.
- Integrate every element/orbital channel over that manifold.
- Divide each element channel by its atom count and apply the absolute
  per-atom thresholds `s=0.4`, `p=1.2`, `d=2.0`, and `f=2.8`.
- Do not cap the number of selected channels or normalize their weights into
  cross-channel fractions.

The next workflow step reads `outputs/projections.txt`, prepares `proj/`, and writes the selected lines into the `begin projections` block of `WANNIER90_WIN`.

## VASP Projection Setup

Use `scripts/prepare_projection.py` after DOS projection analysis succeeds.

Typical command from the project root:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_projection.py" --run-np <actual-mpi-ranks>
```

The script copies `scf` to `proj`, sets `ICHARG = 11` and `NCORE = 1`, removes active `NPAR`, writes `LWANNIER90 = True` and `LWRITE_MMN_AMN = True`, embeds `WANNIER90_WIN`, computes `NUM_WANN`, and sets `NBANDS`.

The workflow supports both scalar and spinor calculations. `NUM_WANN` is computed from `projections.txt` and POSCAR atom counts using `s/p/d/f = 1/3/5/7`; it is doubled only when `LSORBIT` or `LNONCOLLINEAR` is active in the SCF INCAR.
For spinor calculations, projection and band preparation remove the active
`ISPIN` line while preserving `LSORBIT`, `LNONCOLLINEAR`, and three-component
`MAGMOM`; VASP derives the spinor mode from those tags and VASP 6.5 rejects the
legacy `ISPIN = 2` noncollinear combination.

For both `band` and `projection`, ask the user for the actual VASP MPI rank
count and pass it as `--run-np`. Set `NBANDS` from SCF `OUTCAR` as
`2.0 * NELECT`, rounded up to that rank multiple, while also ensuring
`NBANDS >= NUM_WANN` for projection.

## Wannierization Setup

Use `scripts/prepare_wannier.py` after projection succeeds.

Typical command from the project root:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_wannier.py"
```

The script creates a clean `wann` directory and copies only the Wannier90 inputs from `proj`: `wannier90.amn`, `wannier90.eig`, `wannier90.mmn`, and `wannier90.win`. It reads the Fermi energy from `proj/OUTCAR`, reads `proj/EIGENVAL` and `proj/vasprun.xml` to construct the outer disentanglement window, reads `proj/EIGENVAL` again to construct the frozen-window upper limit, caps the frozen window by material type, converts the VASP line-mode `band/KPOINTS` into a Wannier90 `begin kpoint_path` block, writes `bands_num_points = 20`, and updates `wannier90.win`.

The outer disentanglement window is selected from a Fermi-centered ranked
distance:

1. For every k point in `EIGENVAL`, sort band energies by absolute distance to
   `E_F`.
2. Select rank `num_wann + outer_band_margin` at every k point. The default
   `outer_band_margin` is `2`.
3. Take the maximum selected distance across all k points so the window
   contains at least that ranked state everywhere.
4. Write the symmetric window
   `dis_win_min = E_F - max_distance` and
   `dis_win_max = E_F + max_distance`. A `1e-6 eV` numerical padding is
   applied outside the ranked distance so boundary states remain included
   after text serialization.

The frozen-window upper-limit rule is selected from the band energies inside the current outer window:

1. For every k point in `EIGENVAL`, collect all band energies inside the outer window.
2. Sort their absolute distances to `E_F`.
3. Select the `num_wann - froz_band_margin` ranked distance at each k point.
4. Use the smallest selected distance across all k points as the single-sided frozen-window upper-limit half-width.
5. Compare that upper-limit half-width with the default frozen-window cap: `2 eV` for metals and `4 eV` for semiconductors or insulators.
6. Use the smaller half-width, then write `dis_froz_min = E_F - half_width` and `dis_froz_max = E_F + half_width`.

`froz_band_margin` defaults to `1`. Treat the actual Wannier90 `.wout` result as
authoritative even when the initial ranked window satisfies the energy-count
constraint. If `.wout` reports too many states in the frozen window, increase
`--froz-band-margin` by exactly one, regenerate the window, and rerun
Wannier90. Repeat one rank at a time until that error disappears. If `.wout`
reports too few states in the outer window, increase `--outer-band-margin` by
exactly one and repeat in the same way. Read `docs/troubleshooting.md` for the
retry procedure. For large conventional cells that genuinely need a wider
frozen window, override the default cap explicitly with
`--froz-half-width-cap-ev`.

The script writes the production convergence controls:

```text
dis_num_iter = 3000
num_iter = 3000
dis_conv_tol = 1.0e-10
num_cg_steps = 30
guiding_centres = true
search_shells = 25
```

Before running the standalone Wannier stage, follow
`docs/scheduler_guidance.md`: detect Slurm or PBS, then ask the user for the
Wannier90 executable, queue, requested resources, actual MPI ranks, memory,
wall time, and required environment setup.
It also writes `wann/window_summary.json`, which freezes the four
physics-derived initial window edges and `num_wann` for later optimization and
for a candidate-independent fitting assessment window.

It also writes:

```text
write_hr = True
bands_plot = True
bands_num_points = 20
begin kpoint_path
...
end kpoint_path
```

## Final VASP Band From Wannier K-Points

Use `scripts/prepare_band_from_wannier_kpt.py` after Wannier90 has produced `<seedname>_band.kpt`.

Typical command from the project root:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_band_from_wannier_kpt.py" --wann-dir wann --run-np <actual-mpi-ranks>
```

The script reads `wann/wannier90_band.kpt`, writes `band/KPOINTS` as an explicit reciprocal-coordinate VASP KPOINTS file, and updates `band/INCAR` for a non-self-consistent SOC band run. Create and confirm a scheduler job for this final VASP band calculation before plotting or LESA analysis.

## Fitting Plot

Use `scripts/plot_fitting.py` after Wannier90 has produced `<seedname>_band.dat`.

Typical command from the project root:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/plot_fitting.py" --vasp-dir band --wann-dir wann
```

The script reads SOC/noncollinear VASP bands directly from `band/EIGENVAL`, uses the Wannier path coordinate from `wannier90_band.dat`, overlays DFT and Wannier bands, and writes `<project>_fitting.png`. It does not call VASPKIT task `211`. Always produce this PNG for human inspection; LESA text or JSON alone is not a complete fitting deliverable.

## LESA Fitting Quality

Use `scripts/lesa_fitting.py` after Wannier90 has produced `<seedname>_band.dat`
and the final VASP band job has produced `band/EIGENVAL` from `wannier90_band.kpt`.

Typical command from the project root:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/lesa_fitting.py" --config "$AUTO_MLWFS_SKILL_DIR/configuration/fitting.yaml"
```

This is the default agent-readable fitting judgment and should be paired with
the fitting PNG for human inspection. Do not judge fitting quality from the PNG
alone. LESA builds a local relative eigenvalue-spectrum
anchor from DFT states around E_F, searches sliding Wannier windows at each k
point, and classifies anchors as failed, success, or ambiguous. The global
anchor success ratio must pass `decision.min_anchor_success_ratio` before the
initial frozen-window fitting ratio is computed. Treat the LESA anchor as the
coupling/alignment decision and the initial frozen window as the sole fitting
quality decision region.

Read `docs/lesa_fitting.md` for the rule details. The script writes
`outputs/fitting_lesa_summary.json` and `outputs/fitting_lesa_report.txt`.

## Four-Dimensional Window Optimization

If the physics-derived `wann` result does not pass LESA, use
`scripts/optimize_windows.py` with `configuration/optimization.yaml`. Keep
projections, `num_wann`, `NBANDS`, and all non-window Wannier90 inputs fixed.
The only optimization variables are `dis_win_min`, `dis_win_max`,
`dis_froz_min`, and `dis_froz_max`.

Before every `propose`, ask the user how many candidates to run concurrently.
Recommend `10`, but use the user's value. Generate a batch with that explicit
choice:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" --root . propose --batch-size <user-choice>
```

The Python environment used for `propose` requires `numpy`, `PyYAML`, and
`scikit-optimize`; candidate tasks additionally use the existing LESA
dependencies. `propose` checks the BO dependency before creating work.
Ten candidates is a recommendation, not authorization. The search has no
evaluation limit unless `max_evaluations` is configured.
`propose` creates scheduler-independent trial directories and writes
`optimization/batch-XXX.json` with the user-chosen batch size, trial IDs,
candidate windows, work directories, and scheduler inputs still to confirm;
it does not create or submit a job script. Follow `docs/scheduler_guidance.md`,
ask the user for the missing environment and resource values, and create a
Slurm/PBS job or array that invokes `run-trial` once per candidate with its
trial ID and the user-approved Wannier90 command. Submit only after user
confirmation. When every task finishes, collect the batch:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" --root . collect
```

If the reported status remains active, run `propose` again. The optimizer
persists JSON state under `optimization/`, so the loop is restartable without
pickles.

Before a candidate is submitted, require at least `num_wann` outer-window
states at every projection k point, fewer than `num_wann` frozen-window states
at every k point, containment of the frozen window by the outer window, and
the configured frozen minimum width. The outer window may grow or shrink
around the physical guess. The frozen-window search starts from the largest
safe physical guess and explores inward, but it is not a monotonic shrink
sequence.

The objective is the equal-weight normalized sum of the paper-style smoothed
band RMS discrepancy (`sigma = 0.001 eV`) and mean Wannier spread. DFT and
Wannier energies are both referenced to the same DFT Fermi energy, split into
below/near/above-Fermi buckets, and paired by signed distance rank within each
k point. Use the initial physics frozen window read from
`wann/window_summary.json` as the sole LESA fitting region, never a candidate's
moving frozen window. Require at least `0.95` frozen-window state coverage.
LESA is the final usability gate and may stop the search early after a complete
batch. Inspect `optimization/best.json` and still generate the fitting PNG for
the selected candidate.

Every completed optimization trial also writes its own
`optimization/trials/trial-XXXX/fitting.png`. A second copy is archived as
`optimization/plots/trial-XXXX.png`, so plots from earlier batches are never
overwritten and the fitting evolution can be reviewed after the search.

If BO makes no progress toward a LESA pass, do not blindly append points in
the same bounded space. Check whether the candidate with the best band RMS or
highest LESA anchor-success ratio lies on a search boundary. Start a separately
named phase such as `optimization_phase2/`, reuse that candidate as the new
baseline, expand only the boundary-supported directions, and retain the
earlier phase unchanged. When LESA progress and localization conflict, the
follow-up phase may increase the band/LESA contribution and reduce—but not
necessarily remove—the spread contribution.

`examples/agi/optimization_phase2.yaml` is a historical AgI example, not a
default configuration. Replace its baseline trial path with a completed trial
from the user's calculation before using it.

## Scheduler Interaction

Read `docs/scheduler_guidance.md` before creating or submitting any job. Detect
Slurm or PBS from the current environment and existing successful scripts. Ask
the user for the VASP and Wannier90 locations or launch commands, queue,
requested cores, actual MPI ranks, memory, wall time, account, and environment
setup when those values have not already been supplied.

For VASP, read `INCAR`: use the user-specified noncollinear/SOC executable when
`LSORBIT` or `LNONCOLLINEAR` is true, otherwise use the standard executable.
The preparation scripts report this as `vasp_mode`. Generate a scheduler script
for the detected system, show its resolved resources and command, and submit it
only after explicit user confirmation.

If a `band` or `projection` job has credible memory-failure evidence, stop and
ask whether the user wants to reduce actual MPI ranks, increase requested
memory/resources, or do both. If ranks change, regenerate the inputs with the
new `--run-np` before submitting again.

## VASP Band Setup

Use `scripts/vasp_band.py` to prepare a VASP band calculation from `scf/`. Read `docs/band_workflow.md` when changing the band workflow or handling an existing `band/` directory.

Typical command from the project root:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/vasp_band.py" --run-np <actual-mpi-ranks>
```

The script copies `scf` to `band`, deletes inherited `WAVECAR` and run outputs, runs VASPKIT task `303` for bulk by default, and keeps `KPATH.in`/line-mode `KPOINTS` as the path source for Wannier90. If `band/` already exists, the default output becomes `band_agent/`; if both exist, stop and ask the user. Do not submit this preliminary band setup as the final VASP band calculation. After Wannier90 writes `wannier90_band.kpt`, run `prepare_band_from_wannier_kpt.py` to write the explicit final KPOINTS in the active band directory, then create the scheduler job in conversation with the user.

Keep requested scheduler resources separate from actual MPI ranks. For memory
failures, follow `docs/troubleshooting.md` and `docs/scheduler_guidance.md`.

After a VASP band job finishes, run:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/check_success.py" --stage vasp-band --workdir band
```

Proceed only if the script returns success. See `docs/success_criteria.md`.

## Agent-Readable Stage Results

The band-path, projection, Wannier preparation, and final-band preparation
scripts write `stage_summary.json` in the actual output directory. Read this
file instead of parsing console prose. It records the selected directory,
whether the `_agent` suffix was used, physical settings, detected VASP mode,
scheduler inputs to confirm, and the next action.

For newly created `band`, `proj`, and `wann` stages, use the default collision
policy: if the requested directory exists, write the new result to
`<name>_agent`; if both paths exist, stop and ask the user which path to use or
preserve. Carry the `workdir` reported by `stage_summary.json` into later
commands such as `--band-dir`, `--proj-dir`, and `--wann-dir`.
