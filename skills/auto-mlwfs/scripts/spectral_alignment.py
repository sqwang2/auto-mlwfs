#!/usr/bin/env python3
"""Fermi-referenced, same-side rank matching for DFT and Wannier spectra."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SpectralPair:
    bucket: str
    rank: int
    dft_energy_ev: float
    wannier_energy_ev: float | None


def _ranked_buckets(values: np.ndarray, neutral_tolerance_ev: float) -> dict[str, list[float]]:
    energies = [float(value) for value in np.asarray(values, dtype=float)]
    return {
        "below": sorted((value for value in energies if value < -neutral_tolerance_ev), reverse=True),
        "neutral": sorted(
            (value for value in energies if abs(value) <= neutral_tolerance_ev),
            key=lambda value: (abs(value), value),
        ),
        "above": sorted(value for value in energies if value > neutral_tolerance_ev),
    }


def pair_by_signed_fermi_rank(
    dft_energies_ev: np.ndarray,
    wannier_energies_ev: np.ndarray,
    neutral_tolerance_ev: float,
) -> list[SpectralPair]:
    """Pair closest-below with closest-below, and likewise at/above Fermi.

    Both input spectra must already be referenced to the same DFT Fermi energy.
    A missing Wannier rank is represented by ``wannier_energy_ev=None``.
    """

    if neutral_tolerance_ev < 0.0:
        raise ValueError("neutral_tolerance_ev must be non-negative")
    dft = _ranked_buckets(dft_energies_ev, neutral_tolerance_ev)
    wannier = _ranked_buckets(wannier_energies_ev, neutral_tolerance_ev)
    pairs: list[SpectralPair] = []
    for bucket in ("below", "neutral", "above"):
        for rank, dft_energy in enumerate(dft[bucket]):
            wannier_energy = wannier[bucket][rank] if rank < len(wannier[bucket]) else None
            pairs.append(
                SpectralPair(
                    bucket=bucket,
                    rank=rank,
                    dft_energy_ev=dft_energy,
                    wannier_energy_ev=wannier_energy,
                )
            )
    return pairs
