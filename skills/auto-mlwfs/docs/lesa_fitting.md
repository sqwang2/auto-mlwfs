# LESA Fitting Quality

LESA, the local eigenvalue spectrum anchor, is the default machine-readable
criterion for judging Wannier fitting quality after `wannier90_band.dat` exists.
Agents should read the JSON and text outputs from this workflow instead of
judging the fitting PNG visually.

Run from the project root:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/lesa_fitting.py" --config "$AUTO_MLWFS_SKILL_DIR/configuration/fitting.yaml"
```

The script reads SOC/noncollinear VASP bands directly from `band/EIGENVAL`,
uses `band/OUTCAR` for the Fermi level, and compares against Wannier90
`wannier90_band.dat`. VASP `band/KPOINTS` must first be generated from
`wannier90_band.kpt` with `prepare_band_from_wannier_kpt.py`, then the final
VASP band job must be completed. The script writes:

- `outputs/fitting_lesa_summary.json`
- `outputs/fitting_lesa_report.txt`

## Anchor Rule

At each k point, LESA selects `anchor_size` consecutive DFT eigenvalues around
the configured Fermi level. The DFT anchor is converted to a relative local
spectrum:

```text
S_DFT = E_DFT - E_DFT[0]
```

The Wannier spectrum at the same k point is searched by sliding all consecutive
windows of the same size. Each candidate window is converted to the same
relative form and accepted only if every eigenvalue differs by less than
`lesa.tolerance_ev`.

Anchor status is:

- `anchor_failed`: no Wannier window satisfies the tolerance.
- `anchor_success`: exactly one Wannier window satisfies the tolerance.
- `anchor_ambiguous`: multiple Wannier windows satisfy the tolerance.

Only `anchor_success` points are eligible for target-window fitting when
`lesa.require_unique_anchor` is true. Before any target-window fitting is
performed, the global `anchor_success_ratio` must also meet
`decision.min_anchor_success_ratio`, which defaults to `0.95`.

Use the LESA anchor stage to establish a unique DFT/Wannier spectral coupling.
After that coupling passes globally, judge fitting quality only inside the
initial frozen window.

## Window Fitting

After the global anchor gate passes and a unique LESA anchor is found, DFT
states in the initial frozen window are compared
to the same-k Wannier spectrum by set-to-set nearest-energy matching. This is
still a spectrum comparison, not a band-index comparison.

If `anchor_success_ratio < decision.min_anchor_success_ratio`, the script skips
the window fitting stage and does not report window fitting ratios. In that
case the MLWF is not considered to reproduce the DFT band data well enough for
the downstream window comparison to be meaningful.

For each DFT state in a window, the nearest unused Wannier eigenvalue in the
same window is found. The state is counted as fitted when the absolute energy
difference is below that window's tolerance. Fitting ratio is:

```text
matched_states / target_dft_states
```

The sole fitting-region definition is the initial physics frozen window recorded
in `wann/window_summary.json`. It is fixed for baseline and BO trials and never
follows a candidate's moving frozen window.
The frozen-window fitting ratio must be at least `0.95`.

## Decision

The final JSON field `result.status` is `pass` only when all configured
thresholds are satisfied:

- enough k points have successful LESA anchors;
- failed and ambiguous anchors are not excessive;
- frozen-window fitting ratio is above threshold;
- the LESA energy offset standard deviation is not excessive.

The offset check is only a post-anchor sanity check. It is not used during LESA
matching.

## Output Names

The text report uses these field names:

- `result`: final pass/fail status.
- `anchor_success_ratio`: fraction of k points with exactly one LESA anchor.
- `anchor_failed_ratio`: fraction of k points with no LESA anchor.
- `multi_anchor_ratio`: fraction of k points with multiple possible anchors, present only when the global anchor gate passes.
- `shift_energy`: standard deviation of the LESA energy shift across k points, in eV, present only when the global anchor gate passes.
- `fitting_ratio_initial_frozen_window`: frozen-window fitting ratio, present only when the global anchor gate passes.

When the global anchor gate fails, the text report stops after
`anchor_failed_ratio`. The JSON summary still keeps the full diagnostic fields.
