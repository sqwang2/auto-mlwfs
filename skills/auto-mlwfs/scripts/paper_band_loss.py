#!/usr/bin/env python3
"""Paper-style weighted band loss plus Wannier spread parsing."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import numpy as np

from spectral_alignment import pair_by_signed_fermi_rank


FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def _fermi_filter(energy_ev: float, edge_ev: float, sigma_ev: float) -> float:
    argument = (edge_ev - energy_ev) / sigma_ev
    if argument >= 0.0:
        return 1.0 / (1.0 + math.exp(-min(argument, 700.0)))
    exponential = math.exp(max(argument, -700.0))
    return exponential / (1.0 + exponential)


def paper_window_weight(energy_ev: float, minimum_ev: float, maximum_ev: float, sigma_ev: float) -> float:
    if sigma_ev <= 0.0:
        raise ValueError("sigma_ev must be positive")
    return _fermi_filter(energy_ev, maximum_ev, sigma_ev) * (
        1.0 - _fermi_filter(energy_ev, minimum_ev, sigma_ev)
    )


def weighted_band_discrepancy(
    dft_bands_ev: np.ndarray,
    wannier_bands_ev: np.ndarray,
    assessment_min_ev: float,
    assessment_max_ev: float,
    sigma_ev: float = 0.001,
    neutral_tolerance_ev: float = 0.005,
    unmatched_error_ev: float = 5.0,
) -> dict[str, Any]:
    """Return Eq.-(6)-style weighted RMS after signed-distance rank matching."""

    dft = np.asarray(dft_bands_ev, dtype=float)
    wannier = np.asarray(wannier_bands_ev, dtype=float)
    if dft.ndim != 2 or wannier.ndim != 2:
        raise ValueError("DFT and Wannier bands must both be two-dimensional")
    if dft.shape[0] != wannier.shape[0]:
        raise ValueError(f"k-point mismatch: DFT={dft.shape[0]}, Wannier={wannier.shape[0]}")
    if assessment_min_ev >= assessment_max_ev:
        raise ValueError("assessment_min_ev must be smaller than assessment_max_ev")

    weighted_square_error = 0.0
    weight_sum = 0.0
    matched_count = 0
    unmatched_count = 0
    per_kpoint: list[dict[str, int]] = []

    for k_index in range(dft.shape[0]):
        matched_here = 0
        unmatched_here = 0
        for pair in pair_by_signed_fermi_rank(dft[k_index], wannier[k_index], neutral_tolerance_ev):
            energy = pair.dft_energy_ev
            weight = paper_window_weight(energy, assessment_min_ev, assessment_max_ev, sigma_ev)
            if weight <= 1.0e-14:
                continue
            if pair.wannier_energy_ev is None:
                error = float(unmatched_error_ev)
                unmatched_count += 1
                unmatched_here += 1
            else:
                error = float(pair.wannier_energy_ev - energy)
                matched_count += 1
                matched_here += 1
            weighted_square_error += weight * error * error
            weight_sum += weight
        per_kpoint.append({"k_index": k_index, "matched": matched_here, "unmatched": unmatched_here})

    if weight_sum <= 0.0:
        raise ValueError("No DFT states received non-zero weight in the assessment window")
    return {
        "band_rms_ev": float(math.sqrt(weighted_square_error / weight_sum)),
        "weight_sum": float(weight_sum),
        "matched_count": int(matched_count),
        "unmatched_count": int(unmatched_count),
        "per_kpoint": per_kpoint,
    }


def read_final_spreads(wout_path: Path, num_wann: int) -> np.ndarray:
    """Read the last complete set of individual WF spreads from a .wout file."""

    text = wout_path.read_text(encoding="utf-8", errors="ignore")
    pattern = re.compile(
        rf"WF centre and spread\s+\d+\s+\([^)]*\)\s+({FLOAT_PATTERN})",
        flags=re.IGNORECASE,
    )
    values = [float(value) for value in pattern.findall(text)]
    if len(values) < num_wann:
        raise ValueError(
            f"Found only {len(values)} individual spreads in {wout_path}; expected at least {num_wann}"
        )
    return np.asarray(values[-num_wann:], dtype=float)


def spread_summary(wout_path: Path, num_wann: int) -> dict[str, float]:
    spreads = read_final_spreads(wout_path, num_wann)
    return {
        "spread_sum_a2": float(np.sum(spreads)),
        "spread_mean_a2": float(np.mean(spreads)),
        "spread_max_a2": float(np.max(spreads)),
    }
