#!/usr/bin/env python3
"""Prepare the final VASP band run from Wannier90's sampled band k-points."""

from __future__ import annotations

import argparse
from pathlib import Path

from vasp_band import (
    clean_band_dir,
    copy_scf_to_band,
    nbands_from_nelect,
    read_nelect,
    update_incar,
    vasp_mode,
)
from workflow_io import write_stage_summary


def read_wannier_band_kpt(path: Path) -> list[tuple[float, float, float]]:
    if not path.exists():
        raise FileNotFoundError(f"Wannier band k-point file not found: {path}")
    lines = [line.strip() for line in path.read_text(encoding="utf-8", errors="ignore").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"Wannier band k-point file is empty: {path}")
    try:
        expected = int(lines[0].split()[0])
    except ValueError as exc:
        raise ValueError(f"First line of {path} must be the k-point count") from exc

    kpoints: list[tuple[float, float, float]] = []
    for raw_line in lines[1:]:
        parts = raw_line.split()
        if len(parts) < 3:
            continue
        try:
            kpoints.append((float(parts[0]), float(parts[1]), float(parts[2])))
        except ValueError as exc:
            raise ValueError(f"Malformed k-point row in {path}: {raw_line!r}") from exc
    if len(kpoints) != expected:
        raise ValueError(f"{path} declares {expected} k-points but contains {len(kpoints)} rows")
    return kpoints


def write_vasp_explicit_kpoints(path: Path, kpoints: list[tuple[float, float, float]]) -> None:
    lines = [
        "KPOINTS from wannier90_band.kpt",
        str(len(kpoints)),
        "Reciprocal",
    ]
    lines.extend(f"{kx: .12f} {ky: .12f} {kz: .12f} 1.0" for kx, ky, kz in kpoints)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Calculation root containing scf/, band/, and wann/.")
    parser.add_argument("--scf-dir", default="scf")
    parser.add_argument("--band-dir", default="band")
    parser.add_argument("--wann-dir", default="wann")
    parser.add_argument(
        "--run-np",
        type=int,
        required=True,
        help="MPI ranks that will actually run VASP; used to round NBANDS.",
    )
    parser.add_argument("--seedname", default="wannier90")
    parser.add_argument("--existing", choices=["fail", "backup", "overwrite", "continue"], default="continue")
    parser.add_argument(
        "--keep-old-outputs",
        action="store_true",
        help="Do not delete old VASP output files in band/ before writing KPOINTS.",
    )
    args = parser.parse_args()
    if args.run_np <= 0:
        raise ValueError("--run-np must be positive")

    root = Path(args.root).resolve()
    scf_dir = root / args.scf_dir
    band_dir = root / args.band_dir
    wann_dir = root / args.wann_dir

    if not band_dir.exists() or args.existing != "continue":
        copy_scf_to_band(scf_dir, band_dir, args.existing)
    if not args.keep_old_outputs:
        clean_band_dir(band_dir)

    nelect = read_nelect(scf_dir / "OUTCAR")
    nbands = nbands_from_nelect(nelect, args.run_np)
    update_incar(band_dir / "INCAR", nbands=nbands)

    wannier_kpt = wann_dir / f"{args.seedname}_band.kpt"
    kpoints = read_wannier_band_kpt(wannier_kpt)
    write_vasp_explicit_kpoints(band_dir / "KPOINTS", kpoints)
    mode = vasp_mode(band_dir / "INCAR")
    summary_path = write_stage_summary(
        band_dir,
        {
            "schema_version": 1,
            "status": "prepared",
            "stage": "vasp-band-final",
            "workdir": str(band_dir),
            "wannier_dir": str(wann_dir),
            "wannier_kpt": str(wannier_kpt),
            "vasp_kpoints": str(band_dir / "KPOINTS"),
            "kpoints": len(kpoints),
            "nelect": nelect,
            "nbands": nbands,
            "run_np": args.run_np,
            "vasp_mode": mode,
            "scheduler_required": True,
            "scheduler_inputs_to_confirm": [
                "scheduler",
                f"vasp_{mode}_executable",
                "queue_or_partition",
                "requested_nodes_and_cores",
                "memory",
                "walltime",
                "account_and_environment",
            ],
            "next_action": "ask_for_missing_scheduler_inputs_then_request_submission_confirmation",
        },
    )
    print(f"band_dir: {band_dir}")
    print(f"wannier_kpt: {wannier_kpt}")
    print(f"vasp_kpoints: {band_dir / 'KPOINTS'}")
    print(f"kpoints: {len(kpoints)}")
    print(f"nelect: {nelect:g}")
    print(f"nbands: {nbands}")
    print(f"run_np: {args.run_np}")
    print(f"vasp_mode: {mode}")
    print("scheduler_script: user_or_agent_managed")
    print(f"stage_summary: {summary_path}")


if __name__ == "__main__":
    main()
