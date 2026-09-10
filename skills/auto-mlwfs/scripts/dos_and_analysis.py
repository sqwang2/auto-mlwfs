#!/usr/bin/env python3
"""Extract VASP PDOS, choose useful Wannier projections, and plot them."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from pymatgen.io.vasp.outputs import Vasprun
from workflow_io import resolve_config

try:
    import yaml
except Exception:  # pragma: no cover - cluster environments may vary.
    yaml = None


DEFAULT_CONFIG = {
    "paths": {
        "vasprun": "vasprun.xml",
        "output_dir": "outputs",
        "csv_name": "dos_projected.csv",
        "summary_name": "projection_summary.json",
        "recommendations_name": "projections.txt",
        "plot_name": "DOS.png",
    },
    "classification": {"metal_gap_tolerance_ev": 0.05},
    "analysis": {
        "mode": "dos_manifold",
        "zero_density_tolerance": 0.0,
        "manifold_padding_ev": 0.1,
    },
    "selection": {
        "min_integrated_weight": 1.0e-6,
        "absolute_weight_per_atom": {
            "s": 0.4,
            "p": 1.2,
            "d": 2.0,
            "f": 2.8,
        },
    },
    "plot": {
        "x_label": "Energy (eV)",
        "y_label": "PDOS (States/eV)",
        "legend_loc": "best",
        "legend_fs": 15,
        "label_fs": 15,
        "tick_fs": 15,
        "figsize": [6.0, 4.2],
        "dpi": 300,
        "energy_zero": "fermi",
        "xlim": None,
        "ylim": None,
    },
}


def deep_update(base, override):
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def load_config(path: Path) -> dict:
    config = json.loads(json.dumps(DEFAULT_CONFIG))
    if yaml is None:
        raise RuntimeError("PyYAML is required to read projection.yaml.")
    with path.open("r", encoding="utf-8") as handle:
        return deep_update(config, yaml.safe_load(handle) or {})


def spin_label(spin) -> str:
    return getattr(spin, "name", str(spin)).lower()


def integrate_dos(dos, lo: float, hi: float) -> float:
    energies = np.asarray(dos.energies)
    mask = (energies >= lo) & (energies <= hi)
    if mask.sum() < 2:
        return 0.0
    return float(
        sum(np.trapezoid(np.asarray(arr)[mask], energies[mask]) for arr in dos.densities.values())
    )


def get_band_edges(vr: Vasprun) -> dict:
    gap, cbm, vbm, is_direct = vr.eigenvalue_band_properties
    return {
        "gap_ev": float(gap),
        "vbm_ev": float(vbm),
        "cbm_ev": float(cbm),
        "is_direct_gap": bool(is_direct),
        "efermi_ev": float(vr.efermi),
    }


def merge_segments(segments: list[list[float]]) -> list[list[float]]:
    merged: list[list[float]] = []
    for lo, hi in sorted(segments):
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return merged


def dos_manifold_window(complete_dos, edges: dict, config: dict) -> list[dict]:
    """Select the connected DOS manifold around E_F using absolute DOS zeros.

    Each nonzero-DOS segment is expanded on both sides by half the band gap
    plus the configured padding. This follows the reference autoconstruction
    workflow while remaining robust to a configurable numerical zero.
    """

    energies = np.asarray(complete_dos.energies, dtype=float)
    total_density = np.zeros_like(energies)
    for values in complete_dos.densities.values():
        total_density += np.asarray(values, dtype=float)

    analysis = config["analysis"]
    zero_tolerance = float(analysis.get("zero_density_tolerance", 0.0))
    active = total_density > zero_tolerance
    if not bool(np.any(active)):
        raise ValueError("Total DOS contains no states above the configured zero-density tolerance.")

    active_indices = np.flatnonzero(active)
    split_indices = np.where(np.diff(active_indices) > 1)[0] + 1
    groups = np.split(active_indices, split_indices)
    padding = 0.5 * max(float(edges["gap_ev"]), 0.0) + float(
        analysis.get("manifold_padding_ev", 0.1)
    )
    segments = [
        [
            max(float(energies[0]), float(energies[group[0]]) - padding),
            min(float(energies[-1]), float(energies[group[-1]]) + padding),
        ]
        for group in groups
        if group.size
    ]
    merged = merge_segments(segments)
    efermi = float(edges["efermi_ev"])

    containing = [segment for segment in merged if segment[0] <= efermi <= segment[1]]
    if containing:
        chosen = containing[0]
    else:
        below = [segment for segment in merged if segment[1] < efermi]
        above = [segment for segment in merged if segment[0] > efermi]
        if below and above:
            chosen = [below[-1][0], above[0][1]]
        elif below:
            chosen = below[-1]
        elif above:
            chosen = above[0]
        else:  # pragma: no cover - guarded by the active-DOS check above.
            raise ValueError("Could not select a DOS manifold around the Fermi energy.")

    return [
        {
            "name": "dos_manifold",
            "lo": float(chosen[0]),
            "hi": float(chosen[1]),
            "zero_density_tolerance": zero_tolerance,
            "segment_padding_ev": padding,
        }
    ]


def analysis_windows(complete_dos, edges: dict, config: dict) -> list[dict]:
    mode = str(config.get("analysis", {}).get("mode", "dos_manifold")).lower()
    if mode != "dos_manifold":
        raise ValueError(
            f"Unknown projection analysis mode: {mode!r}; only 'dos_manifold' is supported."
        )
    return dos_manifold_window(complete_dos, edges, config)


def collect_channels(complete_dos, windows: list[dict]) -> list[dict]:
    rows = []
    for element in complete_dos.structure.composition.elements:
        atom_count = float(complete_dos.structure.composition[element])
        spd_dos = complete_dos.get_element_spd_dos(element)
        for orbital, dos in spd_dos.items():
            densities = {spin_label(spin): np.asarray(arr) for spin, arr in dos.densities.items()}
            window_weights = {
                window["name"]: integrate_dos(dos, float(window["lo"]), float(window["hi"]))
                for window in windows
            }
            combined_weight = float(sum(window_weights.values()))
            rows.append(
                {
                    "element": str(element),
                    "orbital": orbital.name,
                    "atom_count": atom_count,
                    "dos": dos,
                    "energies": np.asarray(dos.energies),
                    "densities": densities,
                    "window_weights": window_weights,
                    "combined_weight": combined_weight,
                    "absolute_weight_per_atom": combined_weight / atom_count,
                }
            )
    return rows


def select_projections(channels: list[dict], config: dict) -> list[dict]:
    selection = config["selection"]
    min_weight = float(selection["min_integrated_weight"])
    thresholds = {
        str(orbital).lower(): float(value)
        for orbital, value in selection["absolute_weight_per_atom"].items()
    }

    ranked = [c for c in channels if c["combined_weight"] > min_weight]
    ranked.sort(key=lambda c: c["absolute_weight_per_atom"], reverse=True)
    chosen: list[dict] = []
    for channel in ranked:
        orbital = str(channel["orbital"]).lower()
        if orbital not in thresholds:
            raise ValueError(f"No absolute projection threshold configured for orbital {orbital!r}.")
        threshold = thresholds[orbital]
        if channel["absolute_weight_per_atom"] >= threshold:
            item = dict(channel)
            item["selection_threshold"] = threshold
            item["projection"] = f"{item['element']}:{item['orbital']}"
            chosen.append(item)
    return chosen


def write_csv(path: Path, channels: list[dict], efermi: float) -> None:
    window_names = list(channels[0]["window_weights"].keys()) if channels else []
    spin_names = sorted({name for c in channels for name in c["densities"]})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["energy_ev", "energy_minus_fermi_ev", "element", "orbital", "spin", "density"])
        for channel in channels:
            for spin in spin_names:
                density = channel["densities"].get(spin)
                if density is None:
                    continue
                for energy, value in zip(channel["energies"], density):
                    writer.writerow([f"{energy:.10f}", f"{energy - efermi:.10f}", channel["element"], channel["orbital"], spin, f"{value:.10e}"])

    weights_path = path.with_name(path.stem + "_weights.csv")
    with weights_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "element",
            "orbital",
            "atom_count",
            "combined_weight",
            "absolute_weight_per_atom",
        ] + [f"{name}_weight" for name in window_names]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for channel in sorted(channels, key=lambda c: c["combined_weight"], reverse=True):
            row = {
                "element": channel["element"],
                "orbital": channel["orbital"],
                "atom_count": f"{channel['atom_count']:.10g}",
                "combined_weight": f"{channel['combined_weight']:.10e}",
                "absolute_weight_per_atom": f"{channel['absolute_weight_per_atom']:.10e}",
            }
            for name in window_names:
                row[f"{name}_weight"] = f"{channel['window_weights'][name]:.10e}"
            writer.writerow(row)


def plot_selected(path: Path, selected: list[dict], edges: dict, windows: list[dict], config: dict) -> None:
    if not selected:
        return
    import matplotlib.pyplot as plt

    plot_config = config["plot"]
    zero = edges["efermi_ev"] if plot_config.get("energy_zero", "fermi") == "fermi" else 0.0
    fig, ax = plt.subplots(figsize=tuple(plot_config["figsize"]))

    for channel in selected:
        density = sum(channel["densities"].values())
        ax.plot(channel["energies"] - zero, density, label=channel["projection"])

    xlim = plot_config.get("xlim")
    if xlim is None:
        lo = min(window["lo"] for window in windows) - zero
        hi = max(window["hi"] for window in windows) - zero
        pad = max(0.2, 0.15 * (hi - lo))
        xlim = [lo - pad, hi + pad]
    ax.set_xlim(xlim)
    if plot_config.get("ylim") is not None:
        ax.set_ylim(plot_config["ylim"])
    ax.axvline(0.0, color="0.55", lw=0.8, ls="--")
    ax.set_xlabel(plot_config["x_label"], fontsize=int(plot_config["label_fs"]))
    ax.set_ylabel(plot_config["y_label"], fontsize=int(plot_config["label_fs"]))
    ax.tick_params(axis="both", labelsize=int(plot_config["tick_fs"]))
    ax.legend(loc=plot_config["legend_loc"], frameon=False, fontsize=int(plot_config["legend_fs"]))
    fig.tight_layout()
    fig.savefig(path, dpi=int(plot_config["dpi"]))
    plt.close(fig)


def write_summary(
    path: Path,
    selected: list[dict],
    channels: list[dict],
    edges: dict,
    material_type: str,
    windows: list[dict],
    config: dict,
) -> dict:
    total = sum(c["combined_weight"] for c in channels)
    thresholds = config["selection"]["absolute_weight_per_atom"]
    selected_names = {c["projection"] for c in selected}
    summary = {
        "material_type": material_type,
        "band_edges": edges,
        "analysis_mode": config["analysis"]["mode"],
        "analysis_windows": windows,
        "selection_mode": "absolute_integrated_weight_per_atom",
        "absolute_weight_thresholds": thresholds,
        "projections": [c["projection"] for c in selected],
        "ranked_channels": [
            {
                "projection": f"{c['element']}:{c['orbital']}",
                "combined_weight": c["combined_weight"],
                "atom_count": c["atom_count"],
                "absolute_weight_per_atom": c["absolute_weight_per_atom"],
                "selection_threshold": thresholds.get(c["orbital"]),
                "selected": f"{c['element']}:{c['orbital']}" in selected_names,
                "fraction": (c["combined_weight"] / total if total > 0 else 0.0),
                "window_weights": c["window_weights"],
            }
            for c in sorted(channels, key=lambda c: c["combined_weight"], reverse=True)
        ],
    }
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", default=".", help="VASP calculation directory.")
    parser.add_argument("--config", default=None, help="Path to projection.yaml.")
    parser.add_argument("--vasprun", default=None, help="Override vasprun.xml path.")
    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    config_path = resolve_config(Path.cwd(), args.config, "projection.yaml")
    config = load_config(config_path)
    paths = config["paths"]
    vasprun_path = Path(args.vasprun).resolve() if args.vasprun else workdir / paths["vasprun"]
    output_dir = workdir / paths["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)

    vr = Vasprun(str(vasprun_path), parse_potcar_file=False)
    edges = get_band_edges(vr)
    material_type = (
        "metal"
        if edges["gap_ev"] <= float(config["classification"]["metal_gap_tolerance_ev"])
        else "semiconductor"
    )
    windows = analysis_windows(vr.complete_dos, edges, config)
    channels = collect_channels(vr.complete_dos, windows)
    selected = select_projections(channels, config)
    if not selected:
        raise ValueError(
            "No projection channel passed the configured absolute integrated-weight thresholds."
        )

    csv_path = output_dir / paths["csv_name"]
    write_csv(csv_path, channels, edges["efermi_ev"])
    summary = write_summary(
        output_dir / paths["summary_name"],
        selected,
        channels,
        edges,
        material_type,
        windows,
        config,
    )
    (output_dir / paths["recommendations_name"]).write_text("\n".join(summary["projections"]) + "\n", encoding="utf-8")
    plot_selected(output_dir / paths["plot_name"], selected, edges, windows, config)

    print(json.dumps({
        "material_type": material_type,
        "gap_ev": edges["gap_ev"],
        "vbm_ev": edges["vbm_ev"],
        "cbm_ev": edges["cbm_ev"],
        "efermi_ev": edges["efermi_ev"],
        "projections": summary["projections"],
        "output_dir": str(output_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
