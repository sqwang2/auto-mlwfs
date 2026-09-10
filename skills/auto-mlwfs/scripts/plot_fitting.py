#!/usr/bin/env python3
"""Plot DFT and Wannier band fitting on one figure."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


def read_wannier_band(path: Path, fermi_energy: float) -> tuple[np.ndarray, list[np.ndarray]]:
    text = path.read_text(encoding="utf-8", errors="ignore").strip()
    blocks = [block for block in re.split(r"\n\s*\n", text) if block.strip()]
    if not blocks:
        raise ValueError(f"Wannier band file appears empty or malformed: {path}")

    kpoints = np.array([float(line.split()[0]) for line in blocks[0].splitlines() if line.split()])
    bands = []
    for block in blocks:
        values = [float(line.split()[1]) - fermi_energy for line in block.splitlines() if len(line.split()) >= 2]
        bands.append(np.array(values))
    return kpoints, bands


def read_eigenval_bands(path: Path, fermi_energy: float) -> np.ndarray:
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
            bands.append(float(parts[1]))
            index += 1
        values.append([energy - fermi_energy for energy in bands])

    if len(values) != nkpts:
        raise ValueError(f"Expected {nkpts} k-points but read {len(values)} from {path}")
    return np.array(values, dtype=float)


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


def plot_fitting(
    eigenval_file: Path,
    wannier_file: Path,
    output: Path,
    fermi_energy: float,
    ylim: tuple[float, float],
    figsize: tuple[float, float],
) -> None:
    if not eigenval_file.exists():
        raise FileNotFoundError(f"VASP EIGENVAL file not found: {eigenval_file}")
    if not wannier_file.exists():
        raise FileNotFoundError(f"Wannier band file not found: {wannier_file}")

    fig, ax = plt.subplots(figsize=figsize)

    kpoints, bands = read_wannier_band(wannier_file, fermi_energy)
    vasp_bands = read_eigenval_bands(eigenval_file, fermi_energy)
    if vasp_bands.shape[0] != kpoints.size:
        raise ValueError(
            f"k-point count mismatch: EIGENVAL has {vasp_bands.shape[0]}, "
            f"Wannier has {kpoints.size}. Rebuild band/KPOINTS from wannier90_band.kpt."
        )

    ax.plot(kpoints, vasp_bands, linestyle="-", color="black", linewidth=0.6)
    ax.plot(kpoints, np.array(bands).T, color="red", linestyle="--", linewidth=0.6)

    ax.set_ylim(ylim)
    ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("k-point", fontsize=7)
    ax.set_ylabel("Energy (eV)", fontsize=7)
    ax.tick_params(axis="both", which="major", labelsize=8)

    handles = [
        Line2D([0], [0], color="black", linewidth=1.6, linestyle="-"),
        Line2D([0], [0], color="red", linewidth=1.6, linestyle="--"),
    ]
    ax.legend(handles, ["DFT", "Wannier"], loc="best", frameon=True, fontsize=8)

    fig.tight_layout()
    fig.savefig(output, format="png", dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Project root containing band/ and wann/.")
    parser.add_argument("--vasp-dir", default="band", help="VASP band directory containing EIGENVAL and OUTCAR.")
    parser.add_argument("--wann-dir", default="wann")
    parser.add_argument("--seedname", default="wannier90")
    parser.add_argument("--eigenval-file", default=None, help="Defaults to <vasp-dir>/EIGENVAL.")
    parser.add_argument("--wannier-file", default=None, help="Defaults to <wann-dir>/<seedname>_band.dat.")
    parser.add_argument("--fermi-shift", type=float, default=0.0)
    parser.add_argument("--ylim", nargs=2, type=float, default=(-4.0, 4.0), metavar=("YMIN", "YMAX"))
    parser.add_argument("--figsize", nargs=2, type=float, default=(4.0, 3.5), metavar=("WIDTH", "HEIGHT"))
    parser.add_argument("--output", default=None, help="Defaults to <project-name>_fitting.png.")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    vasp_dir = root / args.vasp_dir
    wann_dir = root / args.wann_dir
    eigenval_file = Path(args.eigenval_file).resolve() if args.eigenval_file else vasp_dir / "EIGENVAL"
    wannier_file = Path(args.wannier_file).resolve() if args.wannier_file else wann_dir / f"{args.seedname}_band.dat"
    output = Path(args.output).resolve() if args.output else root / f"{root.name}_fitting.png"

    fermi_energy = read_fermi_from_outcar(vasp_dir / "OUTCAR") + args.fermi_shift
    plot_fitting(
        eigenval_file=eigenval_file,
        wannier_file=wannier_file,
        output=output,
        fermi_energy=fermi_energy,
        ylim=(args.ylim[0], args.ylim[1]),
        figsize=(args.figsize[0], args.figsize[1]),
    )

    print(f"fermi_energy: {fermi_energy:.10f}")
    print(f"eigenval_file: {eigenval_file}")
    print(f"wannier_file: {wannier_file}")
    print(f"output: {output}")


if __name__ == "__main__":
    main()
