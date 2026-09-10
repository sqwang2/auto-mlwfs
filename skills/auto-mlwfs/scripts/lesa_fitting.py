#!/usr/bin/env python3
"""Evaluate Wannier fitting quality with LESA spectrum anchors."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
from workflow_io import resolve_config

try:
    import yaml
except ImportError as exc:  # pragma: no cover - environment guard
    raise SystemExit("PyYAML is required: pip install pyyaml") from exc


def read_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def read_fermi_from_outcar(path: Path) -> float:
    if not path.exists():
        raise FileNotFoundError(f"OUTCAR not found: {path}")
    text = path.read_text(encoding="utf-8", errors="ignore")
    matches = re.findall(r"E-fermi\s*:\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)", text, flags=re.IGNORECASE)
    if not matches:
        matches = re.findall(r"Fermi energy:\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)", text, flags=re.IGNORECASE)
    if not matches:
        raise ValueError(f"Could not find Fermi energy in {path}")
    return float(matches[-1])


def read_eigenval_bands(path: Path, fermi_energy: float) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"VASP EIGENVAL file not found: {path}")
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    header_index = None
    nkpts = nbands = None
    for index in range(5, min(len(lines), 20)):
        parts = lines[index].split()
        if len(parts) >= 3 and all(part.lstrip("+-").isdigit() for part in parts[:3]):
            _nelect, nkpts, nbands = [int(part) for part in parts[:3]]
            header_index = index
            break
    if header_index is None or nkpts is None or nbands is None:
        raise ValueError(f"Could not read NELECT/NKPTS/NBANDS header from {path}")

    values: list[list[float]] = []
    index = header_index + 1
    while index < len(lines) and len(values) < nkpts:
        if not lines[index].strip():
            index += 1
            continue
        if len(lines[index].split()) < 4:
            index += 1
            continue
        index += 1

        bands: list[float] = []
        for _band in range(nbands):
            if index >= len(lines):
                raise ValueError(f"EIGENVAL ended before all bands were read: {path}")
            parts = lines[index].split()
            if len(parts) < 3:
                raise ValueError(f"Malformed EIGENVAL band row: {lines[index]!r}")
            bands.append(float(parts[1]) - fermi_energy)
            index += 1
        values.append(bands)

    if len(values) != nkpts:
        raise ValueError(f"Expected {nkpts} k-points but read {len(values)} from {path}")
    return np.array(values, dtype=float)


def read_wannier_band(path: Path, fermi_energy: float) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"Wannier band file not found: {path}")
    text = path.read_text(encoding="utf-8", errors="ignore").strip()
    blocks = [block for block in re.split(r"\n\s*\n", text) if block.strip()]
    if not blocks:
        raise ValueError(f"Wannier band file appears empty or malformed: {path}")

    kpoints = np.array([float(line.split()[0]) for line in blocks[0].splitlines() if line.split()])
    bands: list[np.ndarray] = []
    for block in blocks:
        values = [float(line.split()[1]) - fermi_energy for line in block.splitlines() if len(line.split()) >= 2]
        if len(values) != len(kpoints):
            raise ValueError(f"Inconsistent k-point count in Wannier band file: {path}")
        bands.append(np.array(values, dtype=float))
    return kpoints, np.array(bands, dtype=float).T


def match_wannier_kpoints(
    dft_kpoints: np.ndarray,
    wann_kpoints: np.ndarray,
    max_distance: float,
) -> tuple[np.ndarray, np.ndarray]:
    indices: list[int] = []
    distances: list[float] = []
    for k_value in dft_kpoints:
        nearest = int(np.argmin(np.abs(wann_kpoints - k_value)))
        distance = float(abs(wann_kpoints[nearest] - k_value))
        if distance > max_distance:
            raise ValueError(
                f"Could not match DFT k={k_value:.8f} to Wannier k path within {max_distance}; "
                f"nearest distance is {distance:.8f}"
            )
        indices.append(nearest)
        distances.append(distance)
    return np.array(indices, dtype=int), np.array(distances, dtype=float)


def select_anchor_indices(eigenvalues: np.ndarray, anchor_size: int, center_ev: float) -> tuple[int, int]:
    if anchor_size <= 0:
        raise ValueError("anchor_size must be positive")
    if eigenvalues.size < anchor_size:
        raise ValueError("not enough DFT eigenvalues for requested anchor_size")

    center_index = int(np.argmin(np.abs(eigenvalues - center_ev)))
    start = center_index - anchor_size // 2
    start = max(0, min(start, eigenvalues.size - anchor_size))
    return start, start + anchor_size


def relative_spectrum(values: np.ndarray) -> np.ndarray:
    return values - values[0]


def find_lesa_candidates(
    dft_anchor: np.ndarray,
    wannier_values: np.ndarray,
    tolerance_ev: float,
) -> list[dict[str, Any]]:
    anchor_size = dft_anchor.size
    dft_relative = relative_spectrum(dft_anchor)
    candidates: list[dict[str, Any]] = []

    for start in range(0, wannier_values.size - anchor_size + 1):
        window = wannier_values[start : start + anchor_size]
        delta = np.abs(dft_relative - relative_spectrum(window))
        if bool(np.all(delta < tolerance_ev)):
            candidates.append(
                {
                    "start": start,
                    "end": start + anchor_size,
                    "max_relative_error_ev": float(np.max(delta)),
                    "offset_ev": float(dft_anchor[0] - window[0]),
                }
            )
    return candidates


def match_energy_regions_ratio(
    dft_values: np.ndarray,
    wannier_values: np.ndarray,
    regions_ev: list[list[float]],
    tolerance_ev: float,
) -> dict[str, Any]:
    def region_mask(values: np.ndarray) -> np.ndarray:
        mask = np.zeros(values.shape, dtype=bool)
        for emin_ev, emax_ev in regions_ev:
            mask |= (values >= float(emin_ev)) & (values <= float(emax_ev))
        return mask

    dft_window = np.sort(dft_values[region_mask(dft_values)])
    wann_window = list(np.sort(wannier_values[region_mask(wannier_values)]))
    matched_errors: list[float] = []

    for value in dft_window:
        if not wann_window:
            break
        distances = np.abs(np.array(wann_window) - value)
        nearest = int(np.argmin(distances))
        error = float(distances[nearest])
        if error <= tolerance_ev:
            matched_errors.append(error)
            wann_window.pop(nearest)

    target_count = int(dft_window.size)
    matched_count = len(matched_errors)
    return {
        "target_count": target_count,
        "matched_count": matched_count,
        "ratio": float(matched_count / target_count) if target_count else None,
        "mean_error_ev": float(np.mean(matched_errors)) if matched_errors else None,
        "max_error_ev": float(np.max(matched_errors)) if matched_errors else None,
    }


def ratio_ok(value: float | None, threshold: float) -> bool:
    return value is not None and value >= threshold


def summarize_window(records: list[dict[str, Any]], name: str) -> dict[str, Any]:
    total_targets = sum(record[name]["target_count"] for record in records)
    total_matched = sum(record[name]["matched_count"] for record in records)
    ratios = [record[name]["ratio"] for record in records if record[name]["ratio"] is not None]
    errors = [
        record[name]["mean_error_ev"]
        for record in records
        if record[name]["mean_error_ev"] is not None
    ]
    return {
        "target_count": int(total_targets),
        "matched_count": int(total_matched),
        "ratio": float(total_matched / total_targets) if total_targets else None,
        "mean_kpoint_ratio": float(np.mean(ratios)) if ratios else None,
        "mean_error_ev": float(np.mean(errors)) if errors else None,
    }


def window_ratio_label(name: str, config: dict[str, Any]) -> str:
    source = config.get("source") or config.get("resolved_source")
    if source == "initial_frozen_window":
        return f"{name}_ratio_{source}"
    emin = float(config.get("emin_ev"))
    emax = float(config.get("emax_ev"))
    return f"{name}_ratio_{emin:g}_to_{emax:g}_ev"


def initial_frozen_window_edges(summary: dict[str, Any]) -> tuple[float, float]:
    nested = summary.get("physics_initial_window", {})
    minimum_key = "dis_froz_min"
    maximum_key = "dis_froz_max"
    initial_minimum_key = "initial_dis_froz_min"
    initial_maximum_key = "initial_dis_froz_max"
    minimum = nested.get(minimum_key, summary.get(initial_minimum_key, summary.get(minimum_key)))
    maximum = nested.get(maximum_key, summary.get(initial_maximum_key, summary.get(maximum_key)))
    if minimum is None or maximum is None:
        raise ValueError("Initial frozen window is missing from window_summary.json")
    return float(minimum), float(maximum)


def resolve_frozen_window(
    root: Path,
    paths: dict[str, Any],
    windows: dict[str, Any],
    reference_fermi_energy: float,
) -> dict[str, Any]:
    frozen = dict(windows.get("frozen_window", {}))
    source = str(frozen.get("source", "initial_frozen_window")).lower()
    summary_path = root / paths.get("window_summary_file", "wann/window_summary.json")
    if not summary_path.exists():
        raise FileNotFoundError(
            f"Physics-derived LESA windows require the window summary: {summary_path}"
        )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    fermi = float(reference_fermi_energy)

    if source != "initial_frozen_window":
        raise ValueError(f"Unknown frozen-window source: {source!r}")
    frozen_min, frozen_max = initial_frozen_window_edges(summary)
    if not frozen_min < frozen_max:
        raise ValueError(
            f"Initial frozen window is not ordered in {summary_path}: "
            f"{frozen_min} >= {frozen_max}"
        )
    frozen["emin_ev"] = frozen_min - fermi
    frozen["emax_ev"] = frozen_max - fermi
    frozen["regions_ev"] = [[frozen["emin_ev"], frozen["emax_ev"]]]
    frozen["resolved_source"] = source
    frozen["window_summary_file"] = str(summary_path)
    if not float(frozen["emin_ev"]) < float(frozen["emax_ev"]):
        raise ValueError(
            f"Frozen window is not ordered: {frozen['emin_ev']} >= {frozen['emax_ev']}"
        )
    return frozen


def analyze(root: Path, config: dict[str, Any]) -> dict[str, Any]:
    paths = config.get("paths", {})
    lesa = config.get("lesa", {})
    windows = config.get("windows", {})
    decision_config = config.get("decision", {})

    eigenval_file = root / paths.get("eigenval_file", "band/EIGENVAL")
    vasp_outcar = root / paths.get("vasp_outcar", "band/OUTCAR")
    wannier_file = root / paths.get("wannier_file", "wann/wannier90_band.dat")
    fermi_energy = read_fermi_from_outcar(vasp_outcar)
    dft_bands = read_eigenval_bands(eigenval_file, fermi_energy)
    wann_kpoints, wann_bands = read_wannier_band(wannier_file, fermi_energy)
    if dft_bands.shape[0] != wann_bands.shape[0]:
        raise ValueError(
            f"k-point count mismatch: EIGENVAL has {dft_bands.shape[0]}, "
            f"Wannier has {wann_bands.shape[0]}. Run prepare_band_from_wannier_kpt.py "
            "and rerun the VASP band calculation."
        )

    anchor_size = int(lesa.get("anchor_size", 10))
    anchor_tolerance = float(lesa.get("tolerance_ev", 0.05))
    fermi_level = float(lesa.get("fermi_level_ev", 0.0))
    require_unique = bool(lesa.get("require_unique_anchor", True))
    min_anchor_success_ratio = float(decision_config.get("min_anchor_success_ratio", 0.95))

    frozen_config = resolve_frozen_window(
        root,
        paths,
        windows,
        reference_fermi_energy=fermi_energy,
    )
    records: list[dict[str, Any]] = []
    anchor_counts = {"anchor_failed": 0, "anchor_success": 0, "anchor_ambiguous": 0}
    offsets: list[float] = []

    for k_index, dft_values in enumerate(dft_bands):
        wann_values = wann_bands[k_index]
        start, end = select_anchor_indices(dft_values, anchor_size, fermi_level)
        dft_anchor = dft_values[start:end]
        candidates = find_lesa_candidates(dft_anchor, wann_values, anchor_tolerance)

        if len(candidates) == 0:
            status = "anchor_failed"
        elif len(candidates) == 1:
            status = "anchor_success"
        else:
            status = "anchor_ambiguous"
        anchor_counts[status] += 1

        record: dict[str, Any] = {
            "k_index": k_index,
            "k_path": float(wann_kpoints[k_index]),
            "wannier_k_index": k_index,
            "wannier_k": float(wann_kpoints[k_index]),
            "k_distance": 0.0,
            "anchor_status": status,
            "dft_anchor_start": start,
            "dft_anchor_end": end,
            "candidate_count": len(candidates),
            "candidates": candidates[:5],
        }

        records.append(record)

    total_kpoints = len(records)
    success_ratio = anchor_counts["anchor_success"] / total_kpoints if total_kpoints else math.nan
    failed_ratio = anchor_counts["anchor_failed"] / total_kpoints if total_kpoints else math.nan
    ambiguous_ratio = anchor_counts["anchor_ambiguous"] / total_kpoints if total_kpoints else math.nan
    anchor_gate_passed = success_ratio >= min_anchor_success_ratio

    if anchor_gate_passed:
        for record in records:
            k_index = int(record["k_index"])
            dft_values = dft_bands[k_index]
            wann_values = wann_bands[k_index]
            candidates = record["candidates"]
            status = record["anchor_status"]
            can_compare = status == "anchor_success" or (not require_unique and candidates)
            if can_compare:
                offsets.append(float(candidates[0]["offset_ev"]))
                if frozen_config.get("enabled", True):
                    record["frozen_window"] = match_energy_regions_ratio(
                        dft_values,
                        wann_values,
                        frozen_config["regions_ev"],
                        float(frozen_config.get("tolerance_ev", 0.08)),
                    )
    comparable_records = [record for record in records if "frozen_window" in record]

    summary: dict[str, Any] = {
        "inputs": {
            "eigenval_file": str(eigenval_file),
            "vasp_outcar": str(vasp_outcar),
            "wannier_file": str(wannier_file),
            "dft_kpoints": int(dft_bands.shape[0]),
            "wannier_kpoints": int(wann_bands.shape[0]),
            "dft_bands": int(dft_bands.shape[1]),
            "wannier_bands": int(wann_bands.shape[1]),
            "fermi_energy": fermi_energy,
            "same_kpoints": True,
        },
        "parameters": {
            "anchor_size": anchor_size,
            "anchor_tolerance_ev": anchor_tolerance,
            "fermi_level_ev": fermi_level,
            "require_unique_anchor": require_unique,
            "min_anchor_success_ratio": min_anchor_success_ratio,
            "frozen_window": {
                "source": frozen_config["resolved_source"],
                "emin_ev": float(frozen_config["emin_ev"]),
                "emax_ev": float(frozen_config["emax_ev"]),
                "regions_ev": frozen_config["regions_ev"],
                "tolerance_ev": float(frozen_config.get("tolerance_ev", 0.08)),
                "window_summary_file": frozen_config.get("window_summary_file"),
            },
        },
        "anchor": {
            **anchor_counts,
            "total_kpoints": total_kpoints,
            "success_ratio": float(success_ratio),
            "failed_ratio": float(failed_ratio),
            "ambiguous_ratio": float(ambiguous_ratio),
        },
        "post_anchor_comparison": {
            "skipped": not anchor_gate_passed,
            "reason": (
                None
                if anchor_gate_passed
                else f"anchor_success_ratio {success_ratio:.6f} is below required {min_anchor_success_ratio:.6f}"
            ),
            "comparable_kpoints": len(comparable_records),
        },
        "offset": {
            "count": len(offsets),
            "mean_ev": float(np.mean(offsets)) if offsets else None,
            "std_ev": float(np.std(offsets)) if offsets else None,
            "min_ev": float(np.min(offsets)) if offsets else None,
            "max_ev": float(np.max(offsets)) if offsets else None,
        },
        "per_kpoint": records,
    }

    if anchor_gate_passed and frozen_config.get("enabled", True):
        summary["frozen_window"] = summarize_window(comparable_records, "frozen_window")

    frozen_ratio_label = window_ratio_label("fitting", frozen_config)
    checks = {
        "anchor_success_ratio": anchor_gate_passed,
        "anchor_failed_ratio": failed_ratio <= float(decision_config.get("max_anchor_failed_ratio", 0.10)),
        "multi_anchor_ratio": ambiguous_ratio <= float(decision_config.get("max_anchor_ambiguous_ratio", 0.20)),
        "post_anchor_comparison": anchor_gate_passed,
        frozen_ratio_label: ratio_ok(
            summary.get("frozen_window", {}).get("ratio"),
            float(decision_config.get("min_frozen_window_fitting_ratio", 0.95)),
        ),
        "shift_energy": (
            summary["offset"]["std_ev"] is not None
            and summary["offset"]["std_ev"] <= float(decision_config.get("max_offset_std_ev", 0.10))
        ),
    }
    summary["result"] = {
        "status": "pass" if all(checks.values()) else "fail",
        "checks": checks,
        "thresholds": decision_config,
    }
    return summary


def write_outputs(root: Path, config: dict[str, Any], summary: dict[str, Any]) -> tuple[Path, Path]:
    paths = config.get("paths", {})
    output_dir = root / paths.get("output_dir", "outputs")
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / paths.get("summary_name", "fitting_lesa_summary.json")
    report_path = output_dir / paths.get("report_name", "fitting_lesa_report.txt")

    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    frozen_config = summary["parameters"]["frozen_window"]
    report_lines = [
        f"result: {summary['result']['status']}",
        f"anchor_size: {summary['parameters']['anchor_size']}",
        f"anchor_tolerance_ev: {summary['parameters']['anchor_tolerance_ev']}",
        f"anchor_success_ratio: {summary['anchor']['success_ratio']:.6f}",
        f"anchor_failed_ratio: {summary['anchor']['failed_ratio']:.6f}",
    ]
    if not summary["post_anchor_comparison"]["skipped"]:
        report_lines.append(f"multi_anchor_ratio: {summary['anchor']['ambiguous_ratio']:.6f}")
        report_lines.append(f"shift_energy: {summary['offset']['std_ev']}")
    if not summary["post_anchor_comparison"]["skipped"] and "frozen_window" in summary:
        report_lines.append(
            f"{window_ratio_label('fitting', frozen_config)}: {summary['frozen_window']['ratio']}"
        )
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    return summary_path, report_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Project root containing band/ and wann/.")
    parser.add_argument("--config", default=None, help="YAML relative to --root; defaults to bundled fitting.yaml.")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    config_path = resolve_config(root, args.config, "fitting.yaml")

    config = read_config(config_path)
    summary = analyze(root, config)
    summary_path, report_path = write_outputs(root, config, summary)

    print(f"result: {summary['result']['status']}")
    print(f"summary: {summary_path}")
    print(f"report: {report_path}")
    print(f"anchor_size: {summary['parameters']['anchor_size']}")
    print(f"anchor_tolerance_ev: {summary['parameters']['anchor_tolerance_ev']}")
    print(f"anchor_success_ratio: {summary['anchor']['success_ratio']:.6f}")
    print(f"anchor_failed_ratio: {summary['anchor']['failed_ratio']:.6f}")
    if not summary["post_anchor_comparison"]["skipped"]:
        print(f"multi_anchor_ratio: {summary['anchor']['ambiguous_ratio']:.6f}")
        print(f"shift_energy: {summary['offset']['std_ev']}")
    frozen_config = summary["parameters"]["frozen_window"]
    if not summary["post_anchor_comparison"]["skipped"] and "frozen_window" in summary:
        print(f"{window_ratio_label('fitting', frozen_config)}: {summary['frozen_window']['ratio']}")
    return 0 if summary["result"]["status"] == "pass" else 1


if __name__ == "__main__":
    sys.exit(main())
