# Troubleshooting

## Initial Wannier90 Window Count Failures

The initial window construction is intentionally aggressive: start from the
largest physically useful frozen window and the smallest physically useful
outer window so the result remains a strong BO starting point. A window may
satisfy the ranked energy-count construction and still fail when Wannier90
actually runs. Do not replace this initial strategy with a more conservative
window. Treat the error reported in `<seedname>.wout` as authoritative and
adjust only the rank used to construct the failing window.

If `.wout` reports that the frozen window contains too many states, increase
`froz_band_margin` by one. This changes the frozen-window construction rank
from `num_wann - m` to `num_wann - (m + 1)` and shrinks the frozen window.
Regenerate the `wann` directory and rerun Wannier90. If the same error remains,
increase the margin by one again. Continue one rank at a time until the frozen
window error disappears:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_wannier.py" --froz-band-margin 2
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_wannier.py" --froz-band-margin 3
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_wannier.py" --froz-band-margin 4
```

If `.wout` reports that the outer window does not contain enough states,
increase `outer_band_margin` by one. This changes the outer-window construction
rank from `num_wann + m` to `num_wann + (m + 1)` and extends the outer window.
Regenerate the `wann` directory and rerun Wannier90. If the same error remains,
increase the margin by one again. Continue one rank at a time until the outer
window error disappears:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_wannier.py" --outer-band-margin 3
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_wannier.py" --outer-band-margin 4
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_wannier.py" --outer-band-margin 5
```

Handle the two failures independently. Preserve the last successful value of
the other margin while adjusting the failing one. Do not skip several ranks at
once, because the first error-free window is the intended aggressive BO
starting point.

## VASP Band Memory Failure

For `vasp-band`, memory pressure is a common failure mode. Confirm the failure
from the scheduler log, `OUTCAR`, or application output. Then ask the user
whether to reduce the actual MPI ranks, increase the explicit memory/resource
request, or do both. Read `scheduler_guidance.md` before changing the job.

Use this when `band` fails with memory-like symptoms in `OUTCAR`, `.error`, or `vasp-band.out`, such as process kill, allocation failure, MPI abort after large arrays, or segmentation faults during diagonalization.

If the actual MPI rank count changes, rebuild the band inputs so `NBANDS` is
rounded consistently:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/vasp_band.py" --run-np <new-actual-mpi-ranks>
```

Regenerate the Slurm or PBS script with the user-selected resources and submit
only after user confirmation.
## `kmesh_get_bvector: Not enough bvectors found`

Symptom in projection NSCF output:

```text
Calling wannier_setup of wannier90 in library mode
kmesh_get_bvector: Not enough bvectors found
Error on node 0: examine the output/error files for details
```

This means Wannier90 could not identify enough b-vectors for the current structure and k mesh under the current tolerance. The automatic retry ladder is:

```text
kmesh_tol = 0.0001
kmesh_tol = 0.0002
kmesh_tol = 0.0005
kmesh_tol = 0.001
```

Do not increase `kmesh_tol` beyond `0.001` automatically. If `0.001` still fails, stop the workflow and ask the user to inspect the initial structure, lattice precision, symmetry, and k mesh.

## Projection Out Of Memory

Symptom in `.error` or `out`:

```text
Detected 1 oom-kill event(s)
task 0: Out Of Memory
KILLED BY SIGNAL: 9 (Killed)
```

For projection jobs, this may happen after Wannier90 setup succeeds, especially with SOC and large `NBANDS`.
Ask the user whether to keep the requested allocation and reduce the actual MPI
ranks, increase the explicit memory/resource request, or use both changes. The
effect of requesting more cores depends on the cluster's memory accounting, so
do not assume it increases available memory. If MPI ranks change, rerun
`prepare_projection.py --run-np <new-actual-mpi-ranks>` before submission.
