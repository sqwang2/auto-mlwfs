#!/usr/bin/env python3
"""Prepare a Wannier90 run from a completed VASP projection directory."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import time
from pathlib import Path

import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun

from workflow_io import agent_output_dir, write_stage_summary



WANNIER_INPUT_SUFFIXES = ("amn", "eig", "mmn", "win")
VASP_TIMING_FOOTER = "General timing and accounting informations for this job"
FROZEN_HALF_WIDTH_CAPS_EV = {
    "metal": 2.0,
    "semiconductor": 4.0,
}


def prepare_wann_dir(proj_dir: Path, wann_dir: Path, seedname: str, mode: str) -> list[str]:
    if not proj_dir.exists():
        raise FileNotFoundError(f"Projection directory not found: {proj_dir}")
    if wann_dir.exists():
        if mode == "fail":
            raise FileExistsError(f"Wannier directory already exists: {wann_dir}")
        if mode == "backup":
            index = 1
            while True:
                backup = wann_dir.with_name(f"{wann_dir.name}.bak{index}")
                if not backup.exists():
                    wann_dir.rename(backup)
                    break
                index += 1
        elif mode == "overwrite":
            try:
                shutil.rmtree(wann_dir)
            except OSError:
                stale = wann_dir.with_name(f"{wann_dir.name}.stale.{int(time.time())}")
                wann_dir.rename(stale)
    wann_dir.mkdir(parents=True, exist_ok=True)

    copied = []
    for suffix in WANNIER_INPUT_SUFFIXES:
        name = f"{seedname}.{suffix}"
        source = proj_dir / name
        target = wann_dir / name
        if not source.exists() or source.stat().st_size == 0:
            raise RuntimeError(f"Projection output is missing or empty: {source}")
        shutil.copy2(source, target)
        copied.append(name)
    return copied


def read_fermi_from_outcar(outcar_path: Path) -> float:
    text = outcar_path.read_text(encoding="utf-8", errors="ignore")
    patterns = [
        r"Fermi energy:\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)",
        r"E-fermi\s*:\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)",
    ]
    for pattern in patterns:
        matches = re.findall(pattern, text, flags=re.IGNORECASE)
        if matches:
            return float(matches[-1])
    raise ValueError(f"Could not find Fermi energy in {outcar_path}")


def total_dos_density(vr: Vasprun) -> np.ndarray:
    densities = vr.tdos.densities
    if Spin.up in densities:
        return np.asarray(densities[Spin.up], dtype=float)
    return sum(np.asarray(values, dtype=float) for values in densities.values())


def pairwise(items: list[float]) -> list[list[float]]:
    return [[items[i], items[i + 1]] for i in range(0, len(items) - 1, 2)]


def merge_windows(windows: list[list[float]], weights: list[float]) -> tuple[list[list[float]], list[float]]:
    merged_windows: list[list[float]] = []
    merged_weights: list[float] = []
    for window, weight in sorted(zip(windows, weights), key=lambda item: item[0][0]):
        if merged_windows and merged_windows[-1][1] >= window[0]:
            merged_windows[-1][1] = max(merged_windows[-1][1], window[1])
            merged_weights[-1] += weight
        else:
            merged_windows.append([window[0], window[1]])
            merged_weights.append(weight)
    return merged_windows, merged_weights


def outer_window_from_zero_dos(vasprun_path: Path, extra_shift: tuple[float, float] = (0.0, 0.0), zero_tol: float = 0.0) -> dict:
    vr = Vasprun(str(vasprun_path), parse_potcar_file=False)
    energies = np.asarray(vr.tdos.energies, dtype=float)
    dos = total_dos_density(vr)
    idos = np.asarray(vr.idos.densities[Spin.up], dtype=float)
    gap, _cbm, _vbm, _is_direct = vr.eigenvalue_band_properties
    efermi = float(vr.efermi)
    extension = (float(gap) + 0.2) / 2.0

    zero = np.abs(dos) <= zero_tol
    boundaries = []
    iboundaries = []
    for index in range(1, len(energies) - 1):
        if zero[index] and not zero[index - 1]:
            boundaries.append(float(energies[index]))
            iboundaries.append(float(idos[index]))
        if zero[index] and not zero[index + 1]:
            boundaries.append(float(energies[index]))
            iboundaries.append(float(idos[index]))

    if len(boundaries) % 2 != 0:
        if not zero[0]:
            boundaries = [float(energies[0])] + boundaries
            iboundaries = [float(idos[0])] + iboundaries
        else:
            boundaries.append(float(energies[-1]))
            iboundaries.append(float(idos[-1]))

    windows = [
        [start - extension + extra_shift[0], end + extension + extra_shift[1]]
        for start, end in pairwise(boundaries)
    ]
    weights = [end - start for start, end in pairwise(iboundaries)]
    windows, weights = merge_windows(windows, weights)

    chosen = []
    electron_count = 0.0
    for window, weight in zip(windows, weights):
        if window[0] <= efermi <= window[1]:
            chosen = window
            electron_count = weight
            break

    if not chosen:
        chosen = [float(energies[0]), float(energies[-1])]
        lower_weight = 0.0
        upper_weight = 0.0
        for window, weight in zip(windows, weights):
            if window[0] > efermi:
                chosen[1] = window[1]
                upper_weight = weight
                break
        for window, weight in reversed(list(zip(windows, weights))):
            if window[1] < efermi:
                chosen[0] = window[0]
                lower_weight = weight
                break
        electron_count = lower_weight + upper_weight

    return {
        "dis_win_min": float(chosen[0]),
        "dis_win_max": float(chosen[1]),
        "gap_ev": float(gap),
        "efermi_vasprun_ev": efermi,
        "zero_dos_electron_count": float(electron_count),
        "material_type": "metal" if float(gap) <= 0.05 else "semiconductor",
    }


def zero_dos_segments(vasprun_path: Path, zero_tol: float = 0.0) -> dict:
    vr = Vasprun(str(vasprun_path), parse_potcar_file=False)
    energies = np.asarray(vr.tdos.energies, dtype=float)
    dos = total_dos_density(vr)
    gap, _cbm, _vbm, _is_direct = vr.eigenvalue_band_properties
    efermi = float(vr.efermi)

    zero = np.abs(dos) <= zero_tol
    segments: list[tuple[float, float]] = []
    start_index = None
    for index, is_zero in enumerate(zero):
        if is_zero and start_index is None:
            start_index = index
        elif not is_zero and start_index is not None:
            segments.append((float(energies[start_index]), float(energies[index - 1])))
            start_index = None
    if start_index is not None:
        segments.append((float(energies[start_index]), float(energies[-1])))

    return {
        "segments": segments,
        "energy_min": float(energies[0]),
        "energy_max": float(energies[-1]),
        "gap_ev": float(gap),
        "efermi_vasprun_ev": efermi,
        "material_type": "metal" if float(gap) <= 0.05 else "semiconductor",
    }


def read_num_wann(win_path: Path) -> int:
    text = win_path.read_text(encoding="utf-8", errors="ignore")
    matches = re.findall(r"^\s*num_wann\s*=\s*(\d+)\s*$", text, flags=re.IGNORECASE | re.MULTILINE)
    if not matches:
        raise ValueError(f"Could not find num_wann in {win_path}")
    return int(matches[-1])


def eigenval_energies_by_kpoint(eigenval_path: Path) -> list[list[float]]:
    lines = eigenval_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 7:
        raise ValueError(f"EIGENVAL is too short: {eigenval_path}")

    header_index = None
    nkpts = nbands = None
    for index in range(5, min(len(lines), 20)):
        line = lines[index]
        parts = line.split()
        if len(parts) >= 3 and all(part.lstrip("+-").isdigit() for part in parts[:3]):
            _nelect, nkpts, nbands = [int(part) for part in parts[:3]]
            header_index = index
            break
    if header_index is None or nkpts is None or nbands is None:
        raise ValueError(f"Could not read NELECT/NKPTS/NBANDS header from {eigenval_path}")

    kpoint_energies: list[list[float]] = []
    index = header_index + 1
    while index < len(lines) and len(kpoint_energies) < nkpts:
        if not lines[index].strip():
            index += 1
            continue
        kpoint_line = lines[index].split()
        if len(kpoint_line) < 4:
            index += 1
            continue
        index += 1

        energies: list[float] = []
        for _band in range(nbands):
            if index >= len(lines):
                raise ValueError(f"EIGENVAL ended before all bands were read: {eigenval_path}")
            parts = lines[index].split()
            if len(parts) < 3:
                raise ValueError(f"Malformed EIGENVAL band row: {lines[index]!r}")
            numeric = [float(part) for part in parts[1:]]
            if len(numeric) >= 4:
                energies.extend(numeric[:2])
            else:
                energies.append(numeric[0])
            index += 1
        kpoint_energies.append(energies)

    if len(kpoint_energies) != nkpts:
        raise ValueError(f"Expected {nkpts} k-points but read {len(kpoint_energies)} from {eigenval_path}")
    return kpoint_energies


def fermi_centered_energy_block(eigenval_path: Path, fermi_energy: float, state_count: int) -> dict:
    if state_count < 1:
        raise ValueError(f"state_count must be positive, got {state_count}")

    selected: list[float] = []
    min_available_states = None
    kpoint_count = 0
    for kpoint_index, energies in enumerate(eigenval_energies_by_kpoint(eigenval_path), start=1):
        ordered = sorted((float(energy) for energy in energies), key=lambda energy: (abs(energy - fermi_energy), energy))
        min_available_states = len(ordered) if min_available_states is None else min(min_available_states, len(ordered))
        if len(ordered) < state_count:
            raise ValueError(
                f"k-point {kpoint_index} has only {len(ordered)} bands, "
                f"but {state_count} are required for the Fermi-centered block."
            )
        selected.extend(ordered[:state_count])
        kpoint_count += 1

    if not selected:
        raise ValueError(f"No k-points were read from {eigenval_path}")

    return {
        "block_min": float(min(selected)),
        "block_max": float(max(selected)),
        "block_state_count": int(state_count),
        "block_kpoint_count": int(kpoint_count),
        "block_min_available_states": int(min_available_states or 0),
    }


def outer_window_from_eigenval_blocks(
    vasprun_path: Path,
    eigenval_path: Path,
    fermi_energy: float,
    num_wann: int,
    outer_band_margin: int = 2,
    zero_tol: float = 0.0,
) -> dict:
    rank = int(num_wann) + int(outer_band_margin)
    if rank <= int(num_wann):
        raise ValueError(
            f"outer_band_margin must be positive, got {outer_band_margin}"
        )

    ranked_distances = []
    min_available_states = None
    kpoint_count = 0
    for kpoint_index, energies in enumerate(
        eigenval_energies_by_kpoint(eigenval_path),
        start=1,
    ):
        distances = sorted(abs(float(energy) - fermi_energy) for energy in energies)
        min_available_states = (
            len(distances)
            if min_available_states is None
            else min(min_available_states, len(distances))
        )
        if len(distances) < rank:
            raise ValueError(
                f"k-point {kpoint_index} has only {len(distances)} bands, "
                f"but outer rank {rank} is required."
            )
        ranked_distances.append(float(distances[rank - 1]))
        kpoint_count += 1

    if not ranked_distances:
        raise ValueError(f"No k-points were read from {eigenval_path}")

    ranked_distance_max = float(max(ranked_distances))
    numerical_padding = 1.0e-6
    half_width = ranked_distance_max + numerical_padding
    dis_win_min = float(fermi_energy - half_width)
    dis_win_max = float(fermi_energy + half_width)
    zero_info = zero_dos_segments(vasprun_path, zero_tol=zero_tol)
    idos_at = idos_interpolator(vasprun_path)
    window_states = float(idos_at(dis_win_max) - idos_at(dis_win_min))
    return {
        "dis_win_min": dis_win_min,
        "dis_win_max": dis_win_max,
        "initial_dis_win_min": dis_win_min,
        "initial_dis_win_max": dis_win_max,
        "initial_window_states": window_states,
        "final_window_states": window_states,
        "target_window_states": float(rank),
        "outer_shrink_iterations": 0,
        "outer_window_mode": "symmetric_ranked_distance",
        "outer_band_margin": int(outer_band_margin),
        "outer_rank": int(rank),
        "outer_half_width": half_width,
        "outer_numerical_padding_ev": numerical_padding,
        "outer_ranked_distance_min": float(min(ranked_distances)),
        "outer_ranked_distance_max": ranked_distance_max,
        "outer_kpoint_count": int(kpoint_count),
        "outer_min_available_states": int(min_available_states or 0),
        "gap_ev": zero_info["gap_ev"],
        "efermi_vasprun_ev": zero_info["efermi_vasprun_ev"],
        "zero_dos_electron_count": window_states,
        "material_type": zero_info["material_type"],
    }


def frozen_window_limit_from_eigenval(
    eigenval_path: Path,
    fermi_energy: float,
    dis_win_min: float,
    dis_win_max: float,
    num_wann: int,
    froz_band_margin: int = 1,
) -> dict:
    rank = int(num_wann) - int(froz_band_margin)
    if rank < 1:
        raise ValueError(f"froz_band_margin={froz_band_margin} leaves no frozen bands for num_wann={num_wann}")

    half_widths = []
    min_outer_states = None
    for kpoint_index, energies in enumerate(eigenval_energies_by_kpoint(eigenval_path), start=1):
        distances = sorted(
            abs(float(energy) - fermi_energy)
            for energy in energies
            if dis_win_min <= float(energy) <= dis_win_max
        )
        min_outer_states = len(distances) if min_outer_states is None else min(min_outer_states, len(distances))
        if len(distances) < rank:
            raise ValueError(
                f"k-point {kpoint_index} has only {len(distances)} bands inside the outer window, "
                f"but {rank} are required by num_wann={num_wann} and froz_band_margin={froz_band_margin}."
            )
        half_widths.append(float(distances[rank - 1]))

    if not half_widths:
        raise ValueError(f"No k-points were read from {eigenval_path}")

    half_width = min(half_widths)
    return {
        "froz_limit_half_width": float(half_width),
        "froz_band_margin": int(froz_band_margin),
        "froz_rank": int(rank),
        "froz_kpoint_count": len(half_widths),
        "froz_min_outer_window_states": int(min_outer_states or 0),
    }


def capped_frozen_window(limit: dict, fermi_energy: float, material_type: str, cap_override: float | None = None) -> dict:
    if material_type not in FROZEN_HALF_WIDTH_CAPS_EV:
        raise ValueError(f"Unknown material_type for frozen-window cap: {material_type!r}")
    cap = float(cap_override) if cap_override is not None else FROZEN_HALF_WIDTH_CAPS_EV[material_type]
    if cap <= 0.0:
        raise ValueError(f"Frozen-window cap must be positive, got {cap}")
    half_width = min(float(limit["froz_limit_half_width"]), cap)
    return {
        **limit,
        "dis_froz_min": float(fermi_energy - half_width),
        "dis_froz_max": float(fermi_energy + half_width),
        "froz_half_width": float(half_width),
        "froz_half_width_cap": float(cap),
        "froz_uncapped_half_width": float(limit["froz_limit_half_width"]),
    }


def idos_interpolator(vasprun_path: Path):
    vr = Vasprun(str(vasprun_path), parse_potcar_file=False)
    energies = np.asarray(vr.idos.energies, dtype=float)
    idos = np.asarray(vr.idos.densities[Spin.up], dtype=float)
    return lambda energy: float(np.interp(energy, energies, idos))


def shrink_outer_window(
    vasprun_path: Path,
    dis_win_min: float,
    dis_win_max: float,
    frozen_min: float,
    frozen_max: float,
    fermi_energy: float,
    num_wann: int,
    max_states_factor: float = 1.2,
    step_ev: float = 0.1,
) -> dict:
    idos_at = idos_interpolator(vasprun_path)
    target_states = max_states_factor * num_wann

    left = float(dis_win_min)
    right = float(dis_win_max)
    initial_states = idos_at(right) - idos_at(left)
    iterations = 0

    while idos_at(right) - idos_at(left) > target_states:
        left_room = float(frozen_min) - left
        right_room = right - float(frozen_max)
        if left_room <= 1.0e-8 and right_room <= 1.0e-8:
            break

        left_distance = abs(left - fermi_energy) if left_room > 1.0e-8 else -1.0
        right_distance = abs(right - fermi_energy) if right_room > 1.0e-8 else -1.0

        if right_distance >= left_distance and right_room > 1.0e-8:
            right = max(right - step_ev, float(frozen_max))
        elif left_room > 1.0e-8:
            left = min(left + step_ev, float(frozen_min))
        else:
            break
        iterations += 1

    final_states = idos_at(right) - idos_at(left)
    return {
        "dis_win_min": left,
        "dis_win_max": right,
        "initial_dis_win_min": float(dis_win_min),
        "initial_dis_win_max": float(dis_win_max),
        "initial_window_states": float(initial_states),
        "final_window_states": float(final_states),
        "target_window_states": float(target_states),
        "outer_shrink_iterations": iterations,
    }


def remove_existing_parameters(win_text: str, keys: set[str]) -> str:
    lines = []
    for line in win_text.splitlines():
        stripped = line.strip()
        if "=" in stripped:
            key = stripped.split("=", 1)[0].strip().lower()
            if key in keys:
                continue
        lines.append(line)
    return "\n".join(lines).rstrip() + "\n"


def remove_block(win_text: str, begin: str, end: str) -> str:
    pattern = re.compile(
        rf"^\s*{re.escape(begin)}\s*$.*?^\s*{re.escape(end)}\s*$\n?",
        flags=re.IGNORECASE | re.MULTILINE | re.DOTALL,
    )
    return re.sub(pattern, "", win_text)


def normalize_kpoint_label(label: str) -> str:
    if label.upper() in {"GAMMA", "GAM"}:
        return "G"
    return label


def parse_vasp_line_mode_kpoints(kpoints_path: Path) -> list[tuple[str, list[float], str, list[float]]]:
    if not kpoints_path.exists():
        raise FileNotFoundError(f"Band KPOINTS file not found: {kpoints_path}")
    all_lines = kpoints_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    rows = []
    in_line_mode = False
    for raw_line in all_lines:
        stripped = raw_line.strip()
        if not stripped:
            continue
        if stripped.lower() == "line-mode":
            in_line_mode = True
            continue
        if not in_line_mode:
            continue
        if stripped.lower() in {"reciprocal", "cartesian"}:
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        try:
            coords = [float(parts[0]), float(parts[1]), float(parts[2])]
        except ValueError:
            continue
        label = normalize_kpoint_label(parts[3])
        rows.append((label, coords))
    if len(rows) < 2:
        raise ValueError(f"No line-mode k-point pairs found in {kpoints_path}")
    if len(rows) % 2 != 0:
        raise ValueError(f"Line-mode KPOINTS has an odd number of endpoints: {kpoints_path}")
    segments = []
    for index in range(0, len(rows), 2):
        start_label, start_coords = rows[index]
        end_label, end_coords = rows[index + 1]
        segments.append((start_label, start_coords, end_label, end_coords))
    return segments


def format_kpoint_path(segments: list[tuple[str, list[float], str, list[float]]]) -> list[str]:
    lines = ["begin kpoint_path"]
    for start_label, start_coords, end_label, end_coords in segments:
        start = " ".join(f"{value:.3f}" for value in start_coords)
        end = " ".join(f"{value:.3f}" for value in end_coords)
        lines.append(f"{start_label} {start} {end_label} {end}")
    lines.append("end kpoint_path")
    return lines


def update_win(
    win_path: Path,
    settings: dict,
    kpoint_segments: list[tuple[str, list[float], str, list[float]]],
    bands_num_points: int,
) -> None:
    keys = {
        "dis_win_min",
        "dis_win_max",
        "dis_num_iter",
        "num_iter",
        "dis_conv_tol",
        "num_cg_steps",
        "guiding_centres",
        "search_shells",
        "dis_froz_min",
        "dis_froz_max",
        "fermi_energy",
        "write_hr",
        "bands_plot",
        "bands_num_points",
    }
    text = win_path.read_text(encoding="utf-8", errors="ignore")
    text = remove_block(text, "begin kpoint_path", "end kpoint_path")
    text = remove_existing_parameters(text, keys)
    insert = [
        f"dis_win_min = {settings['dis_win_min']:.10f}",
        f"dis_win_max = {settings['dis_win_max']:.10f}",
        "",
        f"dis_num_iter = {int(settings['dis_num_iter'])}",
        "num_iter = 3000",
        "dis_conv_tol = 1.0e-10",
        "num_cg_steps = 30",
        "guiding_centres = true",
        "search_shells = 25",
        "",
        f"dis_froz_min = {settings['dis_froz_min']:.10f}",
        f"dis_froz_max = {settings['dis_froz_max']:.10f}",
        "",
        f"fermi_energy = {settings['fermi_energy']:.10f}",
        "",
        "write_hr = True",
        "bands_plot = True",
        f"bands_num_points = {bands_num_points}",
        "",
        *format_kpoint_path(kpoint_segments),
        "",
    ]
    marker = "# This part was generated automatically by VASP"
    if marker in text:
        text = text.replace(marker, "\n".join(insert) + marker, 1)
    else:
        text = text.rstrip() + "\n\n" + "\n".join(insert)
    win_path.write_text(text, encoding="utf-8")


def assert_projection_success(proj_dir: Path, seedname: str) -> None:
    outcar = proj_dir / "OUTCAR"
    if VASP_TIMING_FOOTER not in outcar.read_text(encoding="utf-8", errors="ignore"):
        raise RuntimeError("Projection OUTCAR does not contain VASP timing footer.")
    for suffix in ("amn", "mmn", "eig", "win"):
        path = proj_dir / f"{seedname}.{suffix}"
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f"Projection output is missing or empty: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Calculation root containing proj/ and wann/.")
    parser.add_argument("--proj-dir", default="proj")
    parser.add_argument("--band-dir", default="band", help="Band directory containing the line-mode KPOINTS path.")
    parser.add_argument("--wann-dir", default="wann")
    parser.add_argument("--seedname", default="wannier90")
    parser.add_argument(
        "--existing",
        choices=["agent", "fail", "backup", "overwrite", "continue"],
        default="agent",
        help="On collision, default to <wann-dir>_agent; fail if both paths exist.",
    )
    parser.add_argument("--dis-num-iter", type=int, default=3000)
    parser.add_argument("--froz-band-margin", type=int, default=1)
    parser.add_argument("--froz-half-width-cap-ev", type=float, default=None)
    parser.add_argument(
        "--outer-band-margin",
        type=int,
        default=2,
        help="Use rank num_wann + margin at every k point to set the symmetric outer window.",
    )
    parser.add_argument("--shrink-step-ev", type=float, default=0.1)
    parser.add_argument("--zero-dos-tol", type=float, default=0.0)
    parser.add_argument(
        "--bands-num-points",
        type=int,
        default=20,
        help="Wannier90 band divisions along the first k-path section.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.bands_num_points <= 0:
        raise ValueError("--bands-num-points must be positive")

    root = Path(args.root).resolve()
    proj_dir = root / args.proj_dir
    requested_wann_dir = root / args.wann_dir
    wann_dir = agent_output_dir(requested_wann_dir) if args.existing == "agent" else requested_wann_dir
    existing_mode = "fail" if args.existing == "agent" else args.existing
    band_kpoints = root / args.band_dir / "KPOINTS"
    seedname = args.seedname

    assert_projection_success(proj_dir, seedname)
    fermi_energy = read_fermi_from_outcar(proj_dir / "OUTCAR")
    num_wann = read_num_wann(proj_dir / f"{seedname}.win")
    window = outer_window_from_eigenval_blocks(
        vasprun_path=proj_dir / "vasprun.xml",
        eigenval_path=proj_dir / "EIGENVAL",
        fermi_energy=fermi_energy,
        num_wann=num_wann,
        outer_band_margin=args.outer_band_margin,
        zero_tol=args.zero_dos_tol,
    )
    frozen_limit = frozen_window_limit_from_eigenval(
        eigenval_path=proj_dir / "EIGENVAL",
        fermi_energy=fermi_energy,
        dis_win_min=window["dis_win_min"],
        dis_win_max=window["dis_win_max"],
        num_wann=num_wann,
        froz_band_margin=args.froz_band_margin,
    )
    frozen = capped_frozen_window(
        frozen_limit,
        fermi_energy,
        window["material_type"],
        cap_override=args.froz_half_width_cap_ev,
    )
    settings = {
        **window,
        "dis_num_iter": args.dis_num_iter,
        **frozen,
        "initial_dis_froz_min": frozen["dis_froz_min"],
        "initial_dis_froz_max": frozen["dis_froz_max"],
        "fermi_energy": fermi_energy,
        "num_wann": num_wann,
    }
    kpoint_segments = parse_vasp_line_mode_kpoints(band_kpoints)

    if args.dry_run:
        print(f"wann_dir: {wann_dir}")
        for key in (
            "material_type",
            "gap_ev",
            "fermi_energy",
            "dis_win_min",
            "dis_win_max",
            "dis_froz_min",
            "dis_froz_max",
            "froz_half_width",
            "froz_limit_half_width",
            "froz_half_width_cap",
            "froz_uncapped_half_width",
            "froz_band_margin",
            "froz_rank",
            "froz_kpoint_count",
            "froz_min_outer_window_states",
            "dis_num_iter",
            "num_wann",
            "outer_window_mode",
            "outer_band_margin",
            "outer_rank",
            "outer_half_width",
            "outer_numerical_padding_ev",
            "outer_ranked_distance_min",
            "outer_ranked_distance_max",
            "outer_kpoint_count",
            "outer_min_available_states",
            "initial_dis_win_min",
            "initial_dis_win_max",
            "initial_dis_froz_min",
            "initial_dis_froz_max",
            "initial_window_states",
            "final_window_states",
            "target_window_states",
            "outer_shrink_iterations",
        ):
            print(f"{key}: {settings[key]}")
        print(f"bands_num_points: {args.bands_num_points}")
        print(f"kpoint_path_segments: {len(kpoint_segments)}")
        for line in format_kpoint_path(kpoint_segments):
            print(line)
        print("dry_run: true")
        return

    copied = prepare_wann_dir(proj_dir, wann_dir, seedname, existing_mode)
    update_win(wann_dir / f"{seedname}.win", settings, kpoint_segments, args.bands_num_points)
    window_summary = {
        "schema_version": 1,
        "seedname": seedname,
        "num_wann": int(num_wann),
        "fermi_energy": float(fermi_energy),
        "dis_win_min": float(settings["dis_win_min"]),
        "dis_win_max": float(settings["dis_win_max"]),
        "dis_froz_min": float(settings["dis_froz_min"]),
        "dis_froz_max": float(settings["dis_froz_max"]),
        "physics_initial_window": {
            "dis_win_min": float(settings["dis_win_min"]),
            "dis_win_max": float(settings["dis_win_max"]),
            "dis_froz_min": float(settings["dis_froz_min"]),
            "dis_froz_max": float(settings["dis_froz_max"]),
        },
        "selection_metadata": settings,
    }
    (wann_dir / "window_summary.json").write_text(
        json.dumps(window_summary, indent=2) + "\n",
        encoding="utf-8",
    )
    summary_path = write_stage_summary(
        wann_dir,
        {
            "schema_version": 1,
            "status": "prepared",
            "stage": "wannier",
            "requested_workdir": str(requested_wann_dir),
            "workdir": str(wann_dir),
            "used_agent_suffix": wann_dir != requested_wann_dir,
            "projection_dir": str(proj_dir),
            "band_dir": str(root / args.band_dir),
            "seedname": seedname,
            "num_wann": num_wann,
            "fermi_energy": fermi_energy,
            "windows": {
                "dis_win_min": settings["dis_win_min"],
                "dis_win_max": settings["dis_win_max"],
                "dis_froz_min": settings["dis_froz_min"],
                "dis_froz_max": settings["dis_froz_max"],
            },
            "copied": copied,
            "scheduler_required": True,
            "scheduler_inputs_to_confirm": [
                "scheduler",
                "wannier90_executable",
                "queue_or_partition",
                "requested_nodes_and_cores",
                "actual_mpi_ranks",
                "memory",
                "walltime",
                "account_and_environment",
            ],
            "next_action": "ask_for_missing_scheduler_inputs_then_request_submission_confirmation",
        },
    )
    print(f"wann_dir: {wann_dir}")
    for key in (
        "material_type",
        "gap_ev",
        "fermi_energy",
        "dis_win_min",
        "dis_win_max",
        "dis_froz_min",
        "dis_froz_max",
        "froz_half_width",
        "froz_limit_half_width",
        "froz_half_width_cap",
        "froz_uncapped_half_width",
        "froz_band_margin",
        "froz_rank",
        "froz_kpoint_count",
        "froz_min_outer_window_states",
        "dis_num_iter",
        "num_wann",
        "outer_window_mode",
        "outer_band_margin",
        "outer_rank",
        "outer_half_width",
        "outer_numerical_padding_ev",
        "outer_ranked_distance_min",
        "outer_ranked_distance_max",
        "outer_kpoint_count",
        "outer_min_available_states",
        "initial_dis_win_min",
        "initial_dis_win_max",
        "initial_dis_froz_min",
        "initial_dis_froz_max",
        "initial_window_states",
        "final_window_states",
        "target_window_states",
        "outer_shrink_iterations",
    ):
        print(f"{key}: {settings[key]}")
    print(f"copied: {', '.join(copied)}")
    print(f"window_summary: {wann_dir / 'window_summary.json'}")
    print("scheduler_script: user_or_agent_managed")
    print(f"stage_summary: {summary_path}")


if __name__ == "__main__":
    main()
