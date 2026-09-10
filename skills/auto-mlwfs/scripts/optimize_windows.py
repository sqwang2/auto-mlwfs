#!/usr/bin/env python3
"""Batch Bayesian optimization of four Wannier90 energy-window edges."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

try:
    import yaml
except ImportError as exc:  # pragma: no cover - environment guard
    raise SystemExit("PyYAML is required to run window optimization.") from exc

from lesa_fitting import (
    analyze as lesa_analyze,
    read_config as read_lesa_config,
    read_eigenval_bands,
    read_fermi_from_outcar,
    read_wannier_band,
    write_outputs as write_lesa_outputs,
)
from paper_band_loss import spread_summary, weighted_band_discrepancy
from plot_fitting import plot_fitting
from workflow_io import resolve_config


PARAMETERS = ("dis_win_min", "dis_win_max", "dis_froz_min", "dis_froz_max")
WANNIER_DONE_MARKER = "All done: wannier90 exiting"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def configured_max_evaluations(config: dict[str, Any]) -> int | None:
    value = config.get("bo", {}).get("max_evaluations")
    if value is None:
        return None
    limit = int(value)
    return limit if limit > 0 else None


def atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def append_event(path: Path, event: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps({"time": utc_now(), **event}) + "\n")


def render_placeholders(text: str, context: dict[str, Any]) -> str:
    def replace(match: re.Match[str]) -> str:
        key = match.group(1).strip()
        if key not in context:
            raise KeyError(f"Unknown template value {key!r} in {text!r}")
        return str(context[key])

    return re.sub(r"{{\s*([^}]+)\s*}}", replace, text)


def read_win_values(path: Path) -> dict[str, float | int]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    values: dict[str, float | int] = {}
    for key in (*PARAMETERS, "fermi_energy"):
        matches = re.findall(
            rf"^\s*{key}\s*=\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)",
            text,
            flags=re.IGNORECASE | re.MULTILINE,
        )
        if matches:
            values[key] = float(matches[-1])
    matches = re.findall(r"^\s*num_wann\s*=\s*(\d+)", text, flags=re.IGNORECASE | re.MULTILINE)
    if matches:
        values["num_wann"] = int(matches[-1])
    missing = [key for key in (*PARAMETERS, "num_wann") if key not in values]
    if missing:
        raise ValueError(f"Missing {', '.join(missing)} in {path}")
    return values


def update_win_values(path: Path, candidate: dict[str, float]) -> None:
    text = path.read_text(encoding="utf-8", errors="ignore")
    for key in PARAMETERS:
        pattern = re.compile(
            rf"^(\s*{key}\s*=\s*)([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)(.*)$",
            flags=re.IGNORECASE | re.MULTILINE,
        )
        replacement = rf"\g<1>{candidate[key]:.10f}\g<3>"
        text, count = pattern.subn(replacement, text, count=1)
        if count != 1:
            raise ValueError(f"Could not update {key} in {path}")
    path.write_text(text, encoding="utf-8")


def read_projection_energies(eigenval_path: Path) -> list[list[float]]:
    """Read the projection EIGENVAL using the same spin handling as preparation."""

    lines = eigenval_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    header_index = None
    nkpts = nbands = None
    for index in range(5, min(len(lines), 20)):
        parts = lines[index].split()
        if len(parts) >= 3 and all(part.lstrip("+-").isdigit() for part in parts[:3]):
            _nelect, nkpts, nbands = [int(part) for part in parts[:3]]
            header_index = index
            break
    if header_index is None or nkpts is None or nbands is None:
        raise ValueError(f"Could not read EIGENVAL header from {eigenval_path}")

    result: list[list[float]] = []
    index = header_index + 1
    while index < len(lines) and len(result) < nkpts:
        if not lines[index].strip() or len(lines[index].split()) < 4:
            index += 1
            continue
        index += 1
        energies: list[float] = []
        for _ in range(nbands):
            parts = lines[index].split()
            numeric = [float(value) for value in parts[1:]]
            energies.extend(numeric[:2] if len(numeric) >= 4 else numeric[:1])
            index += 1
        result.append(energies)
    if len(result) != nkpts:
        raise ValueError(f"Expected {nkpts} k-points but read {len(result)} from {eigenval_path}")
    return result


def optimization_paths(root: Path, config: dict[str, Any]) -> dict[str, Path]:
    paths = config.get("paths", {})
    return {
        "projection": root / paths.get("projection_dir", "proj"),
        "band": root / paths.get("band_dir", "band"),
        "baseline": root / paths.get("baseline_wannier_dir", "wann"),
        "optimization": root / paths.get("optimization_dir", "optimization"),
        "fitting_config": resolve_config(root, paths.get("fitting_config"), "fitting.yaml"),
    }


def baseline_metadata(root: Path, config: dict[str, Any]) -> dict[str, Any]:
    paths = optimization_paths(root, config)
    seedname = str(config.get("seedname", "wannier90"))
    summary_path = paths["baseline"] / "window_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    else:
        summary = read_win_values(paths["baseline"] / f"{seedname}.win")
    if "fermi_energy" not in summary:
        summary["fermi_energy"] = read_fermi_from_outcar(paths["projection"] / "OUTCAR")
    return summary


def make_bounds(metadata: dict[str, Any], config: dict[str, Any]) -> dict[str, list[float]]:
    search = config.get("search", {})
    base = {key: float(metadata[key]) for key in PARAMETERS}
    initial_width = base["dis_froz_max"] - base["dis_froz_min"]
    minimum_width = min(
        initial_width,
        max(
            float(search.get("min_frozen_width_ev", 0.2)),
            float(search.get("min_frozen_width_fraction", 0.5)) * initial_width,
        ),
    )
    frozen_min_upper = min(
        base["dis_froz_min"] + float(search.get("frozen_min_inward_ev", 1.0)),
        base["dis_froz_max"] - minimum_width,
    )
    frozen_max_lower = max(
        base["dis_froz_max"] - float(search.get("frozen_max_inward_ev", 1.0)),
        base["dis_froz_min"] + minimum_width,
    )
    return {
        "dis_win_min": [
            base["dis_win_min"] - float(search.get("outer_min_radius_ev", 1.5)),
            base["dis_win_min"] + float(search.get("outer_min_radius_ev", 1.5)),
        ],
        "dis_win_max": [
            base["dis_win_max"] - float(search.get("outer_max_radius_ev", 1.5)),
            base["dis_win_max"] + float(search.get("outer_max_radius_ev", 1.5)),
        ],
        "dis_froz_min": [base["dis_froz_min"], frozen_min_upper],
        "dis_froz_max": [frozen_max_lower, base["dis_froz_max"]],
    }


def candidate_key(candidate: dict[str, float], tolerance: float) -> tuple[int, ...]:
    return tuple(int(round(float(candidate[name]) / tolerance)) for name in PARAMETERS)


def validate_candidate(
    candidate: dict[str, float],
    metadata: dict[str, Any],
    projection_energies: list[list[float]],
    config: dict[str, Any],
) -> tuple[bool, list[str], dict[str, int]]:
    reasons: list[str] = []
    wmin, wmax = candidate["dis_win_min"], candidate["dis_win_max"]
    fmin, fmax = candidate["dis_froz_min"], candidate["dis_froz_max"]
    fermi = float(metadata["fermi_energy"])
    num_wann = int(metadata["num_wann"])
    search = config.get("search", {})
    constraints = config.get("constraints", {})
    initial_width = float(metadata["dis_froz_max"]) - float(metadata["dis_froz_min"])
    minimum_width = min(
        initial_width,
        max(
            float(search.get("min_frozen_width_ev", 0.2)),
            float(search.get("min_frozen_width_fraction", 0.5)) * initial_width,
        ),
    )

    if not wmin < wmax:
        reasons.append("outer window is not ordered")
    if not fmin < fmax:
        reasons.append("frozen window is not ordered")
    if wmin > fmin or wmax < fmax:
        reasons.append("outer window does not contain frozen window")
    if fmax - fmin < minimum_width - 1.0e-12:
        reasons.append(f"frozen width is below {minimum_width:.6f} eV")
    if constraints.get("require_fermi_inside_frozen", True) and not fmin <= fermi <= fmax:
        reasons.append("DFT Fermi energy is outside frozen window")

    outer_counts = [sum(wmin <= energy <= wmax for energy in row) for row in projection_energies]
    frozen_counts = [sum(fmin <= energy <= fmax for energy in row) for row in projection_energies]
    min_outer = min(outer_counts) if outer_counts else 0
    max_frozen = max(frozen_counts) if frozen_counts else 0
    if min_outer < num_wann:
        reasons.append(f"outer window has only {min_outer} states at its worst k point; need >= {num_wann}")
    if max_frozen >= num_wann:
        reasons.append(f"frozen window has {max_frozen} states at its worst k point; need < {num_wann}")
    return not reasons, reasons, {
        "min_outer_states_per_k": int(min_outer),
        "max_frozen_states_per_k": int(max_frozen),
    }


def physical_seed_candidates(metadata: dict[str, Any], bounds: dict[str, list[float]]) -> list[dict[str, float]]:
    base = {name: float(metadata[name]) for name in PARAMETERS}
    seeds = [base]
    changes = [
        {"dis_win_min": 0.5 * (bounds["dis_win_min"][0] + base["dis_win_min"])},
        {"dis_win_max": 0.5 * (bounds["dis_win_max"][1] + base["dis_win_max"])},
        {
            "dis_win_min": 0.5 * (bounds["dis_win_min"][1] + base["dis_win_min"]),
            "dis_win_max": 0.5 * (bounds["dis_win_max"][0] + base["dis_win_max"]),
        },
        {"dis_froz_min": 0.5 * (bounds["dis_froz_min"][1] + base["dis_froz_min"])},
        {"dis_froz_max": 0.5 * (bounds["dis_froz_max"][0] + base["dis_froz_max"])},
        {
            "dis_froz_min": 0.5 * (bounds["dis_froz_min"][1] + base["dis_froz_min"]),
            "dis_froz_max": 0.5 * (bounds["dis_froz_max"][0] + base["dis_froz_max"]),
        },
    ]
    for change in changes:
        seeds.append({**base, **change})
    return seeds


def trial_raw_metrics(
    root: Path,
    trial_dir: Path,
    metadata: dict[str, Any],
    config: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    paths = optimization_paths(root, config)
    seedname = str(config.get("seedname", "wannier90"))
    fermi = float(metadata["fermi_energy"])
    dft = read_eigenval_bands(paths["band"] / "EIGENVAL", fermi)
    _kpath, wannier = read_wannier_band(trial_dir / f"{seedname}_band.dat", fermi)
    loss_config = config.get("loss", {})
    assessment_min = float(metadata["dis_froz_min"]) - fermi
    assessment_max = float(metadata["dis_froz_max"]) - fermi
    discrepancy = weighted_band_discrepancy(
        dft,
        wannier,
        assessment_min,
        assessment_max,
        sigma_ev=float(loss_config.get("sigma_ev", 0.001)),
        neutral_tolerance_ev=(
            float(loss_config.get("sigma_ev", 0.001))
            * float(loss_config.get("neutral_bucket_sigma_multiplier", 5.0))
        ),
        unmatched_error_ev=float(loss_config.get("unmatched_error_ev", 5.0)),
    )
    raw = {
        **discrepancy,
        **spread_summary(trial_dir / f"{seedname}.wout", int(metadata["num_wann"])),
        "assessment_window_ev_relative_to_fermi": [assessment_min, assessment_max],
    }

    plot_path = trial_dir / "fitting.png"
    plot_fitting(
        eigenval_file=paths["band"] / "EIGENVAL",
        wannier_file=trial_dir / f"{seedname}_band.dat",
        output=plot_path,
        fermi_energy=fermi,
        ylim=(-4.0, 4.0),
        figsize=(4.0, 3.5),
    )
    plots_dir = paths["optimization"] / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    archived_plot = plots_dir / f"{trial_dir.name}.png"
    shutil.copy2(plot_path, archived_plot)
    raw["fitting_plot"] = str(plot_path.relative_to(root))
    raw["archived_fitting_plot"] = str(archived_plot.relative_to(root))

    lesa_summary = None
    if config.get("lesa", {}).get("enabled", True):
        lesa_config = read_lesa_config(paths["fitting_config"])
        relative_trial = trial_dir.relative_to(root)
        lesa_config.setdefault("paths", {})["wannier_file"] = str(relative_trial / f"{seedname}_band.dat")
        lesa_config["paths"]["output_dir"] = str(relative_trial / "lesa")
        lesa_summary = lesa_analyze(root, lesa_config)
        write_lesa_outputs(root, lesa_config, lesa_summary)
        raw["lesa_anchor_success_ratio"] = float(
            lesa_summary.get("anchor", {}).get("success_ratio", 0.0)
        )
    return raw, lesa_summary


def normalize_loss(raw: dict[str, Any], baseline_raw: dict[str, Any], config: dict[str, Any]) -> dict[str, float]:
    loss_config = config.get("loss", {})
    band_reference = max(
        float(baseline_raw["band_rms_ev"]),
        float(loss_config.get("normalization_floor_band_ev", 0.001)),
    )
    spread_reference = max(
        float(baseline_raw["spread_mean_a2"]),
        float(loss_config.get("normalization_floor_spread_a2", 0.001)),
    )
    band_normalized = float(raw["band_rms_ev"]) / band_reference
    spread_normalized = float(raw["spread_mean_a2"]) / spread_reference
    anchor_failure = 1.0 - float(raw.get("lesa_anchor_success_ratio", 0.0))
    anchor_failure_weight = float(loss_config.get("lesa_anchor_failure_weight", 0.0))
    total = (
        float(loss_config.get("band_weight", 0.5)) * band_normalized
        + float(loss_config.get("spread_weight", 0.5)) * spread_normalized
        + anchor_failure_weight * anchor_failure
    )
    return {
        "band_normalized": float(band_normalized),
        "spread_normalized": float(spread_normalized),
        "lesa_anchor_failure": float(anchor_failure),
        "lesa_anchor_failure_weight": anchor_failure_weight,
        "total": float(total),
    }


def ensure_state(root: Path, config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Path]]:
    paths = optimization_paths(root, config)
    opt_dir = paths["optimization"]
    state_path = opt_dir / "state.json"
    if state_path.exists():
        return json.loads(state_path.read_text(encoding="utf-8")), paths

    metadata = baseline_metadata(root, config)
    bounds = make_bounds(metadata, config)
    state: dict[str, Any] = {
        "schema_version": 1,
        "created_at": utc_now(),
        "updated_at": utc_now(),
        "status": "active",
        "baseline": metadata,
        "bounds": bounds,
        "baseline_raw": None,
        "trials": [],
        "next_trial_number": 0,
        "next_batch_number": 0,
        "seed_queue": physical_seed_candidates(metadata, bounds),
    }
    opt_dir.mkdir(parents=True, exist_ok=True)

    seedname = str(config.get("seedname", "wannier90"))
    baseline_outputs = (
        paths["baseline"] / f"{seedname}.wout",
        paths["baseline"] / f"{seedname}_band.dat",
    )
    baseline_complete = all(path.exists() and path.stat().st_size > 0 for path in baseline_outputs)
    if baseline_complete:
        baseline_complete = WANNIER_DONE_MARKER in baseline_outputs[0].read_text(
            encoding="utf-8", errors="ignore"
        )
    if baseline_complete:
        raw, lesa = trial_raw_metrics(root, paths["baseline"], metadata, config)
        state["baseline_raw"] = raw
        loss = normalize_loss(raw, raw, config)
        state["trials"].append(
            {
                "trial_id": "baseline-reused",
                "batch": -1,
                "origin": "baseline_reused",
                "candidate": {name: float(metadata[name]) for name in PARAMETERS},
                "status": "completed",
                "raw": raw,
                "loss": loss,
                "lesa_status": lesa["result"]["status"] if lesa else None,
                "completed_at": utc_now(),
            }
        )
        state["seed_queue"] = state["seed_queue"][1:]
        if (
            lesa
            and lesa["result"]["status"] == "pass"
            and config.get("lesa", {}).get("stop_on_pass", True)
        ):
            state["status"] = "stopped"
            state["stop_reason"] = "reused physical baseline already passes LESA"
    atomic_write_json(state_path, state)
    append_event(opt_dir / "trials.jsonl", {"event": "initialized", "baseline_reused": state["baseline_raw"] is not None})
    return state, paths


def save_state(state: dict[str, Any], paths: dict[str, Path]) -> None:
    state["updated_at"] = utc_now()
    atomic_write_json(paths["optimization"] / "state.json", state)


def prepare_trial_directory(
    root: Path,
    paths: dict[str, Path],
    trial_id: str,
    candidate: dict[str, float],
    config: dict[str, Any],
) -> Path:
    trial_dir = paths["optimization"] / "trials" / trial_id
    trial_dir.mkdir(parents=True, exist_ok=False)
    seedname = str(config.get("seedname", "wannier90"))
    mode = str(config.get("execution", {}).get("input_mode", "symlink"))
    for suffix in ("amn", "eig", "mmn"):
        source = paths["baseline"] / f"{seedname}.{suffix}"
        target = trial_dir / source.name
        if not source.exists() or source.stat().st_size == 0:
            raise FileNotFoundError(f"Missing baseline Wannier input: {source}")
        if mode == "symlink":
            try:
                target.symlink_to(os.path.relpath(source, trial_dir))
            except OSError:
                shutil.copy2(source, target)
        else:
            shutil.copy2(source, target)
    shutil.copy2(paths["baseline"] / f"{seedname}.win", trial_dir / f"{seedname}.win")
    update_win_values(trial_dir / f"{seedname}.win", candidate)
    atomic_write_json(trial_dir / "candidate.json", candidate)
    return trial_dir


def completed_observations(state: dict[str, Any], failed_penalty: float) -> tuple[list[list[float]], list[float]]:
    points: list[list[float]] = []
    losses: list[float] = []
    for trial in state["trials"]:
        if trial["status"] not in {"completed", "failed"}:
            continue
        points.append([float(trial["candidate"][name]) for name in PARAMETERS])
        losses.append(float(trial.get("loss", {}).get("total", failed_penalty)))
    return points, losses


def bo_candidates(
    state: dict[str, Any],
    count: int,
    config: dict[str, Any],
    projection_energies: list[list[float]],
) -> list[dict[str, float]]:
    try:
        from skopt import Optimizer
        from skopt.space import Real
    except ImportError as exc:
        raise RuntimeError(
            "scikit-optimize is required after the physics seed queue is exhausted. "
            "Load/install a Python environment containing skopt."
        ) from exc

    bo = config.get("bo", {})
    loss_config = config.get("loss", {})
    dimensions = [
        Real(state["bounds"][name][0], state["bounds"][name][1], name=name)
        for name in PARAMETERS
    ]
    optimizer = Optimizer(
        dimensions=dimensions,
        base_estimator="GP",
        acq_func=str(bo.get("acquisition_function", "EI")),
        n_initial_points=0,
        random_state=int(bo.get("random_seed", 20260724)),
    )
    points, losses = completed_observations(
        state, float(loss_config.get("failed_trial_penalty", 100.0))
    )
    if not points:
        raise RuntimeError("BO requires at least one completed physics seed observation")
    optimizer.tell(points, losses)

    tolerance = float(config.get("search", {}).get("coordinate_tolerance_ev", 0.001))
    existing = {
        candidate_key(trial["candidate"], tolerance)
        for trial in state["trials"]
    }
    proposals: list[dict[str, float]] = []
    pool_factor = int(config.get("search", {}).get("proposal_pool_factor", 12))
    attempts = 0
    while len(proposals) < count and attempts < max(5, count * 4):
        need = count - len(proposals)
        pool = optimizer.ask(
            n_points=max(need, need * pool_factor),
            strategy=str(bo.get("batch_strategy", "cl_min")),
        )
        for vector in pool:
            candidate = {name: float(value) for name, value in zip(PARAMETERS, vector)}
            key = candidate_key(candidate, tolerance)
            valid, _reasons, _counts = validate_candidate(
                candidate, state["baseline"], projection_energies, config
            )
            if valid and key not in existing:
                proposals.append(candidate)
                existing.add(key)
                if len(proposals) == count:
                    break
            else:
                optimizer.tell(vector, float(loss_config.get("failed_trial_penalty", 100.0)))
        attempts += 1
    if len(proposals) < count:
        raise RuntimeError(f"Only found {len(proposals)} valid unique BO candidates; requested {count}")
    return proposals


def require_bo_dependency() -> None:
    try:
        import skopt  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "scikit-optimize is required for this workflow. Install/import skopt "
            "in the Python environment used to propose candidates."
        ) from exc


def propose(root: Path, config: dict[str, Any], batch_size: int) -> int:
    require_bo_dependency()
    if batch_size <= 0:
        raise ValueError("--batch-size must be positive")
    state, paths = ensure_state(root, config)
    if state["status"] != "active":
        print(f"optimization_status: {state['status']}")
        print(f"stop_reason: {state.get('stop_reason')}")
        return 0
    pending = [trial for trial in state["trials"] if trial["status"] == "pending"]
    if pending:
        raise RuntimeError("A proposed batch is still pending; run it and collect results before proposing again")

    max_evaluations = configured_max_evaluations(config)
    count = int(batch_size)
    if max_evaluations is not None:
        remaining = max_evaluations - len(state["trials"])
        if remaining <= 0:
            state["status"] = "stopped"
            state["stop_reason"] = "maximum evaluations reached"
            save_state(state, paths)
            return 0
        count = min(count, remaining)
    projection_energies = read_projection_energies(paths["projection"] / "EIGENVAL")
    tolerance = float(config.get("search", {}).get("coordinate_tolerance_ev", 0.001))
    existing = {candidate_key(trial["candidate"], tolerance) for trial in state["trials"]}
    proposals: list[tuple[dict[str, float], str]] = []

    while state["seed_queue"] and len(proposals) < count:
        candidate = state["seed_queue"].pop(0)
        key = candidate_key(candidate, tolerance)
        valid, _reasons, _counts = validate_candidate(
            candidate, state["baseline"], projection_energies, config
        )
        if valid and key not in existing:
            proposals.append((candidate, "physics_seed"))
            existing.add(key)
    if len(proposals) < count:
        proposals.extend(
            (candidate, "bayesian_optimization")
            for candidate in bo_candidates(
                state, count - len(proposals), config, projection_energies
            )
        )

    batch = int(state["next_batch_number"])
    trial_ids: list[str] = []
    candidate_records: list[dict[str, Any]] = []
    for candidate, origin in proposals:
        trial_id = f"trial-{int(state['next_trial_number']):04d}"
        valid, reasons, counts = validate_candidate(
            candidate, state["baseline"], projection_energies, config
        )
        if not valid:
            raise RuntimeError(f"Internal error: invalid proposal {candidate}: {reasons}")
        prepare_trial_directory(root, paths, trial_id, candidate, config)
        trial = {
            "trial_id": trial_id,
            "batch": batch,
            "origin": origin,
            "candidate": candidate,
            "constraint_counts": counts,
            "status": "pending",
            "proposed_at": utc_now(),
        }
        state["trials"].append(trial)
        append_event(paths["optimization"] / "trials.jsonl", {"event": "proposed", **trial})
        trial_ids.append(trial_id)
        candidate_records.append(
            {
                "trial_id": trial_id,
                "workdir": str(paths["optimization"] / "trials" / trial_id),
                "origin": origin,
                "candidate": candidate,
                "constraint_counts": counts,
            }
        )
        state["next_trial_number"] += 1

    state["next_batch_number"] += 1
    save_state(state, paths)
    batch_summary_path = paths["optimization"] / f"batch-{batch:03d}.json"
    atomic_write_json(
        batch_summary_path,
        {
            "schema_version": 1,
            "status": "proposed",
            "stage": "wannier-window-optimization",
            "batch": batch,
            "requested_batch_size": batch_size,
            "candidate_count": len(trial_ids),
            "trial_ids": trial_ids,
            "candidates": candidate_records,
            "scheduler_required": True,
            "scheduler_inputs_to_confirm": [
                "scheduler",
                "wannier90_executable",
                "queue_or_partition",
                "requested_nodes_and_cores_per_trial",
                "actual_mpi_ranks_per_trial",
                "memory_per_trial",
                "walltime",
                "account_and_environment",
            ],
            "next_action": "create_scheduler_jobs_then_request_submission_confirmation",
        },
    )
    print(f"batch: {batch}")
    print(f"candidate_count: {len(trial_ids)}")
    print(f"trial_ids: {', '.join(trial_ids)}")
    print(f"batch_summary: {batch_summary_path}")
    print("scheduler_script: user_or_agent_managed")
    print("next: create a Slurm/PBS job or array that invokes run-trial once per trial_id")
    return 0


def run_trial(root: Path, config: dict[str, Any], trial_id: str, wannier_command: str) -> int:
    state, paths = ensure_state(root, config)
    trial = next((item for item in state["trials"] if item["trial_id"] == trial_id), None)
    if trial is None:
        raise KeyError(f"Unknown trial id: {trial_id}")
    trial_dir = paths["optimization"] / "trials" / trial_id
    result_path = trial_dir / "result.json"
    if result_path.exists():
        print(f"result_exists: {result_path}")
        return 0

    context = {"seedname": config.get("seedname", "wannier90")}
    command = render_placeholders(wannier_command, context)
    stdout_name = render_placeholders(
        str(config.get("execution", {}).get("wannier_stdout", "{{ seedname }}.wout")),
        context,
    )
    result: dict[str, Any] = {
        "trial_id": trial_id,
        "candidate": trial["candidate"],
        "started_at": utc_now(),
        "command": command,
    }
    try:
        with (trial_dir / stdout_name).open("w", encoding="utf-8") as output:
            completed = subprocess.run(
                command,
                cwd=trial_dir,
                shell=True,
                stdout=output,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode != 0:
            raise RuntimeError(f"Wannier90 exited with code {completed.returncode}")
        wout_text = (trial_dir / stdout_name).read_text(encoding="utf-8", errors="ignore")
        if WANNIER_DONE_MARKER not in wout_text:
            raise RuntimeError(f"Wannier90 output lacks completion marker: {WANNIER_DONE_MARKER}")
        raw, lesa = trial_raw_metrics(root, trial_dir, state["baseline"], config)
        result.update(
            {
                "status": "completed",
                "raw": raw,
                "lesa_status": lesa["result"]["status"] if lesa else None,
            }
        )
    except Exception as exc:
        result.update({"status": "failed", "error": f"{type(exc).__name__}: {exc}"})
    result["completed_at"] = utc_now()
    atomic_write_json(result_path, result)
    print(f"status: {result['status']}")
    print(f"result: {result_path}")
    return 0 if result["status"] == "completed" else 1


def collect(root: Path, config: dict[str, Any]) -> int:
    state, paths = ensure_state(root, config)
    changed = 0
    newly_collected: set[str] = set()
    failed_penalty = float(config.get("loss", {}).get("failed_trial_penalty", 100.0))
    for trial in state["trials"]:
        if trial["status"] != "pending":
            continue
        result_path = paths["optimization"] / "trials" / trial["trial_id"] / "result.json"
        if not result_path.exists():
            continue
        result = json.loads(result_path.read_text(encoding="utf-8"))
        newly_collected.add(trial["trial_id"])
        trial["status"] = result["status"]
        trial["completed_at"] = result.get("completed_at", utc_now())
        trial["lesa_status"] = result.get("lesa_status")
        if result["status"] == "completed":
            trial["raw"] = result["raw"]
            if state["baseline_raw"] is None:
                baseline_candidate = {name: float(state["baseline"][name]) for name in PARAMETERS}
                tolerance = float(config.get("search", {}).get("coordinate_tolerance_ev", 0.001))
                if candidate_key(trial["candidate"], tolerance) == candidate_key(baseline_candidate, tolerance):
                    state["baseline_raw"] = result["raw"]
        else:
            trial["error"] = result.get("error")
            trial["loss"] = {"total": failed_penalty}
        changed += 1

    if any(trial["status"] == "pending" for trial in state["trials"]):
        save_state(state, paths)
        print(f"collected: {changed}")
        print("status: waiting_for_batch")
        return 2
    if state["baseline_raw"] is None:
        successful = [trial for trial in state["trials"] if trial["status"] == "completed"]
        if not successful:
            raise RuntimeError("No trial completed successfully; cannot normalize the loss")
        reference = min(successful, key=lambda trial: (int(trial["batch"]), trial["trial_id"]))
        state["baseline_raw"] = reference["raw"]
        state["normalization_reference_trial_id"] = reference["trial_id"]

    for trial in state["trials"]:
        if trial["trial_id"] not in newly_collected:
            continue
        if trial["status"] == "completed":
            trial["loss"] = normalize_loss(trial["raw"], state["baseline_raw"], config)
            append_event(
                paths["optimization"] / "trials.jsonl",
                {
                    "event": "collected",
                    "trial_id": trial["trial_id"],
                    "status": trial["status"],
                    "loss": trial["loss"],
                    "lesa_status": trial.get("lesa_status"),
                },
            )
        elif trial["status"] == "failed":
            append_event(
                paths["optimization"] / "trials.jsonl",
                {
                    "event": "collected",
                    "trial_id": trial["trial_id"],
                    "status": "failed",
                    "loss": trial["loss"],
                    "error": trial.get("error"),
                },
            )

    completed = [trial for trial in state["trials"] if trial["status"] == "completed"]
    lesa_passed = [trial for trial in completed if trial.get("lesa_status") == "pass"]
    selection_pool = lesa_passed or completed
    if selection_pool:
        if lesa_passed:
            selection = "lesa_passed"
            best = min(
                selection_pool,
                key=lambda trial: (
                    float(trial["loss"]["total"]),
                    float(trial["raw"]["band_rms_ev"]),
                ),
            )
        else:
            selection = "lowest_loss"
            best = min(
                selection_pool,
                key=lambda trial: (
                    float(trial["loss"]["total"]),
                    float(trial["candidate"]["dis_win_max"] - trial["candidate"]["dis_win_min"]),
                    float(trial["raw"]["spread_mean_a2"]),
                ),
            )
        best_record = {
            "selection": selection,
            **best,
        }
        atomic_write_json(paths["optimization"] / "best.json", best_record)
        state["best_trial_id"] = best["trial_id"]

    max_evaluations = configured_max_evaluations(config)
    if lesa_passed and config.get("lesa", {}).get("stop_on_pass", True):
        state["status"] = "stopped"
        state["stop_reason"] = "at least one candidate passed LESA"
    elif max_evaluations is not None and len(state["trials"]) >= max_evaluations:
        state["status"] = "stopped"
        state["stop_reason"] = "maximum evaluations reached"
    save_state(state, paths)
    print(f"collected: {changed}")
    print(f"optimization_status: {state['status']}")
    print(f"best_trial: {state.get('best_trial_id')}")
    if state["status"] == "active":
        print("next: run propose with the same --root and --config and an explicit --batch-size")
    else:
        print(f"stop_reason: {state.get('stop_reason')}")
    return 0


def show_status(root: Path, config: dict[str, Any]) -> int:
    state, paths = ensure_state(root, config)
    counts: dict[str, int] = {}
    for trial in state["trials"]:
        counts[trial["status"]] = counts.get(trial["status"], 0) + 1
    print(f"optimization_status: {state['status']}")
    print(f"trial_counts: {json.dumps(counts, sort_keys=True)}")
    print(f"best_trial: {state.get('best_trial_id')}")
    print(f"state: {paths['optimization'] / 'state.json'}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Calculation root containing proj/, band/, and wann/.")
    parser.add_argument("--config", default=None, help="YAML relative to --root; defaults to bundled optimization.yaml.")
    subparsers = parser.add_subparsers(dest="action", required=True)
    propose_parser = subparsers.add_parser(
        "propose",
        help="Initialize if needed and create the next scheduler-independent candidate batch.",
    )
    propose_parser.add_argument(
        "--batch-size",
        type=int,
        required=True,
        help="Number of candidates explicitly chosen by the user; recommend 10.",
    )
    run_parser = subparsers.add_parser("run-trial", help="Run one candidate inside a user- or agent-managed scheduler job.")
    run_parser.add_argument("--trial-id", required=True)
    run_parser.add_argument(
        "--wannier-command",
        required=True,
        help="Command chosen for the current environment; may contain {{ seedname }}.",
    )
    subparsers.add_parser("collect", help="Collect a finished batch, normalize loss, and update BO history.")
    subparsers.add_parser("status", help="Show persistent optimization status.")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    config_path = resolve_config(root, args.config, "optimization.yaml")
    config = load_yaml(config_path)
    try:
        config["_config_cli_path"] = str(config_path.relative_to(root))
    except ValueError:
        config["_config_cli_path"] = str(config_path)

    if args.action == "propose":
        return propose(root, config, args.batch_size)
    if args.action == "run-trial":
        return run_trial(root, config, args.trial_id, args.wannier_command)
    if args.action == "collect":
        return collect(root, config)
    return show_status(root, config)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"ERROR {type(exc).__name__}: {exc}", file=sys.stderr)
        sys.exit(2)
