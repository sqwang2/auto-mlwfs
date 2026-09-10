# VASP Band Workflow

Use this workflow when the user asks to prepare a SOC/noncollinear VASP band-structure calculation from an SCF result for the auto-MLWF flow.

## Default Procedure

1. Start from a project root containing `scf/`; locate the installed skill separately as described in `SKILL.md`.
2. Prepare `band/` by copying `scf/`.
3. Remove inherited run outputs from `band/`; delete `WAVECAR` by default.
4. Ask the user for the MPI ranks that will actually run VASP, then pass that
   value with `--run-np`. Edit `band/INCAR` minimally:
   - search `ICHARG` and set it to `11`; add it if missing.
   - add `NEDOS = 3001` only if `NEDOS` is missing.
   - set `NBANDS` from SCF `OUTCAR` as `2.0 * NELECT`, rounded up to the `run_np` multiple.
5. Run VASPKIT in `band/` to generate k-path:
   - bulk/3D: task `303`
   - 2D: task `302`
   - 1D: task `301`
6. Keep `KPATH.in` and line-mode `KPOINTS` as the high-symmetry path source for Wannier90.
7. Do not run the preliminary line-mode band job as the final fitting input.
8. After Wannier90 finishes, run `prepare_band_from_wannier_kpt.py` to write an explicit `band/KPOINTS` from `wannier90_band.kpt`.
9. Follow `scheduler_guidance.md` to create the final scheduler script, and
   submit it only after user confirmation.

## Final Explicit KPOINTS

Wannier90 controls the sampled band path through `bands_num_points` and writes the actual points to `<seedname>_band.kpt`. Use those points as the single source of truth for the final VASP band run:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_band_from_wannier_kpt.py" --wann-dir wann --run-np <actual-mpi-ranks>
```

The generated VASP file has this shape:

```text
KPOINTS from wannier90_band.kpt
<N>
Reciprocal
kx ky kz 1.0
...
```

This avoids VASPKIT postprocessing mismatches. Fitting plots and LESA read `band/EIGENVAL` directly.

## Band Resource Decisions

Keep scheduler resources separate from actual MPI ranks:

- requested nodes/cores reserve resources through the current scheduler;
- `--run-np` tells the preparation script how many MPI ranks will actually run
  VASP and controls `NBANDS` rounding;
- memory requests use the syntax and accounting rules of the detected cluster.

For memory-heavy band jobs, ask the user whether to reduce actual MPI ranks or
increase the memory/resource request. See `scheduler_guidance.md` and
`troubleshooting.md`.

## Existing Band Directory

By default, if `band/` exists, the new output is written to `band_agent/`. If
both paths exist, stop and ask the user which directory to use or preserve.
Read the actual `workdir` from `stage_summary.json` and pass it to downstream
commands. Explicit alternatives remain available:

- `--existing backup`: rename old `band/` to `band.bak1`, then rebuild from `scf/`.
- `--existing continue`: reuse current `band/` and only clean/update it.
- `--existing overwrite`: explicitly remove and rebuild the requested directory.

## VASP Binary

Resolve the VASP binary from `band/INCAR`:

- Use `vasp_command_ncl` if `LSORBIT` or `LNONCOLLINEAR` is true.
- Otherwise use `vasp_command_std`.

Ask the user for the corresponding executable paths or launch commands before
creating the scheduler script. The preparation script prints the detected
`vasp_mode`.

Do not use a gamma-only VASP command for this workflow.
