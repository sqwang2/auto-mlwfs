# Scheduler and Resource Guidance

Use this guidance whenever a VASP, Wannier90, or BO trial must run through a
batch scheduler. The skill does not ship a site-specific submit configuration
or generate scheduler scripts. The agent adapts each script to the current
environment in conversation with the user.

## Detect the Scheduler

Before asking for scheduler-specific values, inspect the environment:

- Slurm is indicated by commands such as `sbatch`, `squeue`, and `scontrol`.
- PBS is indicated by commands such as `qsub`, `qstat`, and `pbsnodes`.
- Also inspect successful job scripts already present in the project.
- If both schedulers appear available, neither is available, or the evidence is
  ambiguous, ask the user which scheduler to use.

Do not assume that a scheduler detected on a login host is configured the same
way on another cluster.

## Ask for Site-Specific Inputs

Ask the user for missing values before generating the first job script:

- VASP standard and noncollinear/SOC executable paths or launch commands;
- Wannier90 executable path or launch command;
- queue or partition;
- requested nodes and CPU cores;
- MPI ranks that should actually run the application;
- memory and wall-time request;
- account/project and any required module or environment setup.

Reuse values the user already supplied in the current workflow. Do not ask
again unless the stage requires a different value or an observed failure makes
the earlier choice questionable.

## Select the VASP Executable

Read the active `LSORBIT` and `LNONCOLLINEAR` values from the stage `INCAR`:

- if either is true, use the user-specified noncollinear/SOC executable;
- otherwise use the user-specified standard executable.

The preparation scripts print `vasp_mode: ncl` or `vasp_mode: std`. Treat this
as the preferred machine-readable result. Ask the user if `INCAR` is missing or
ambiguous.

## Keep Allocation and Run Ranks Distinct

The requested CPU allocation and actual MPI ranks are separate decisions.
Pass the actual VASP MPI rank count to input-preparation commands with
`--run-np`; the scripts use it when rounding `NBANDS`.

For example, a job may request more scheduler resources than the number of MPI
ranks used by VASP when that is how the current cluster provides sufficient
memory. Do not assume this improves memory availability: first establish
whether memory is allocated per node, per job, or per CPU on the current site.

## Respond to Memory Failures Interactively

When a `band` or `projection` job shows credible out-of-memory symptoms, stop
and present the evidence to the user. Ask the user to choose between:

1. keeping the requested allocation and reducing the actual MPI ranks; or
2. increasing the explicit memory request or resource allocation.

The user may choose both. Regenerate the job script and, if actual VASP ranks
changed, rerun the relevant preparation command with the new `--run-np` so
`NBANDS` remains consistent. Do not resubmit automatically.

Typical evidence includes scheduler OOM events, signal 9, allocation failures,
or MPI aborts following large allocations. A generic segmentation fault is not
by itself proof of insufficient memory; inspect the surrounding logs.

## BO Jobs

Before `optimize_windows.py propose`, ask the user how many candidates the
batch should contain. Recommend `10`, but use the user's value with
`propose --batch-size <value>`. The command creates scheduler-independent
candidate directories and writes `optimization/batch-XXX.json`. Read its
candidate list and trial work directories instead of parsing console text.
After obtaining the scheduler,
Wannier90 command, queue, resources, and environment setup from the user,
create either individual jobs or a scheduler array. Each task must invoke:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" --root . run-trial \
  --trial-id <trial-id> \
  --wannier-command '<launcher> <wannier90.x> {{ seedname }}'
```

Use the Python executable available in the selected runtime environment. After
all tasks in the batch finish, run `collect`.

## Submission Boundary

Generating and checking a scheduler script is allowed. Display the resolved
queue, allocation, run ranks, memory, wall time, executable, and launch command
before submission. Submit only after the user explicitly confirms that job.
Record the submitted script and job ID with the calculation outputs.
