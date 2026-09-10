# Success Criteria

Use `scripts/check_success.py` as the workflow gate before moving to the next automated step.

## VASP Band

A `vasp-band` calculation is successful only if all hard checks pass:

1. `OUTCAR` contains:

```text
General timing and accounting informations for this job
```

2. `INCAR` contains an active:

```text
ICHARG = 11
```

3. `KPOINTS` is either the preliminary line path:

```text
Line-Mode
```

or the explicit reciprocal-coordinate list generated from
`wannier90_band.kpt` for the final fitting calculation.

Run from the project root:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/check_success.py" --stage vasp-band --workdir band
```

Exit codes:

- `0`: all checks passed.
- `1`: one or more success checks failed.
- `2`: bad usage or missing work directory.

Do not proceed to Wannier fitting if this check fails.

## Projection

A `projection` calculation is successful only if all hard checks pass:

1. `OUTCAR` contains:

```text
General timing and accounting informations for this job
```

2. `<seedname>.wout`, normally `wannier90.wout`, contains:

```text
Finished setting up k-point neighbours.
Exiting wannier_setup in wannier90
```

3. The following files exist and are non-empty:

```text
wannier90.amn
wannier90.mmn
wannier90.eig
```

Run from the project root:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/check_success.py" --stage projection --workdir proj
```

Use `--seedname` if the Wannier90 seed is not `wannier90`.

Do not use `.error` as a hard success or failure criterion for this stage, because it can contain stale output from earlier failed runs. Use `.error`, `out`, and scheduler logs only for troubleshooting.

## Wannier90

A `wannier` calculation is successful only if `<seedname>.wout` contains:

```text
All done: wannier90 exiting
```

and `<seedname>_band.dat` exists and is non-empty.

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/check_success.py" --stage wannier --workdir wann
```
