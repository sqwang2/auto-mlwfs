#!/usr/bin/env python3
"""Prepare a VASP-Wannier90 projection calculation from SCF outputs."""

from __future__ import annotations

import argparse
import re
import shutil
import time
from pathlib import Path

from vasp_band import (
    incar_flags,
    nbands_from_nelect,
    read_nelect,
    vasp_mode,
)
from workflow_io import agent_output_dir, write_stage_summary


ORBITAL_DIMS = {"s": 1, "p": 3, "d": 5, "f": 7}


def copy_scf_to_proj(scf_dir: Path, proj_dir: Path, mode: str) -> None:
    if not scf_dir.exists():
        raise FileNotFoundError(f"SCF directory not found: {scf_dir}")
    if proj_dir.exists():
        if mode == "fail":
            raise FileExistsError(f"Projection directory already exists: {proj_dir}")
        if mode == "backup":
            index = 1
            while True:
                backup = proj_dir.with_name(f"{proj_dir.name}.bak{index}")
                if not backup.exists():
                    proj_dir.rename(backup)
                    break
                index += 1
        elif mode == "overwrite":
            try:
                shutil.rmtree(proj_dir)
            except OSError:
                stale = proj_dir.with_name(f"{proj_dir.name}.stale.{int(time.time())}")
                proj_dir.rename(stale)
        elif mode == "continue":
            return
    shutil.copytree(scf_dir, proj_dir)


def parse_poscar_counts(poscar_path: Path) -> dict[str, int]:
    lines = poscar_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 7:
        raise ValueError(f"POSCAR is too short: {poscar_path}")
    symbols = lines[5].split()
    count_line_index = 6
    if all(re.fullmatch(r"[+-]?\d+", item) for item in symbols):
        raise ValueError("POSCAR must use VASP 5 element symbols for automatic NUM_WANN.")
    counts = [int(item) for item in lines[count_line_index].split()]
    if len(symbols) != len(counts):
        raise ValueError(f"POSCAR element/count mismatch: {symbols} vs {counts}")
    return dict(zip(symbols, counts))


def read_projections(path: Path) -> list[str]:
    projections = []
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        projections.append(line)
    if not projections:
        raise ValueError(f"No projections found in {path}")
    return projections


def compute_num_wann(projections: list[str], atom_counts: dict[str, int], soc: bool = True) -> int:
    total = 0
    for projection in projections:
        if ":" not in projection:
            raise ValueError(f"Projection must look like 'Fe:d': {projection}")
        element, orbital = [part.strip() for part in projection.split(":", 1)]
        orbital_key = orbital.lower()
        if element not in atom_counts:
            raise ValueError(f"Projection element {element!r} is not present in POSCAR.")
        if orbital_key not in ORBITAL_DIMS:
            raise ValueError(f"Unsupported projection orbital {orbital!r}; expected s/p/d/f.")
        total += atom_counts[element] * ORBITAL_DIMS[orbital_key]
    return total * (2 if soc else 1)


def active_incar_value(incar_path: Path, key: str) -> str | None:
    target = key.upper()
    for raw_line in incar_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.split("#", 1)[0].split("!", 1)[0].strip()
        if not line or "=" not in line:
            continue
        lhs, rhs = line.split("=", 1)
        if lhs.strip().upper() == target:
            return rhs.strip().split()[0] if rhs.strip() else ""
    return None


def incar_bool(value: str | None) -> bool:
    return value is not None and value.strip().upper() in {"T", ".TRUE.", "TRUE", "1", "YES"}


def strip_wannier90_win(lines: list[str], start_index: int) -> int:
    index = start_index + 1
    while index < len(lines):
        if '"' in lines[index]:
            return index + 1
        index += 1
    return index


def make_wannier_win(projections: list[str], kmesh_tol: float) -> list[str]:
    lines = ['WANNIER90_WIN       = "begin projections']
    for projection in projections:
        element, orbital = [part.strip() for part in projection.split(":", 1)]
        lines.append(f"{element}: {orbital}")
    lines.append("end projections")
    lines.append("")
    lines.append(f"kmesh_tol = {kmesh_tol:g}")
    lines.append('"')
    return lines


def update_incar(incar_path: Path, projections: list[str], num_wann: int, nbands: int, kmesh_tol: float) -> list[str]:
    spinor = incar_bool(active_incar_value(incar_path, "LSORBIT")) or incar_bool(
        active_incar_value(incar_path, "LNONCOLLINEAR")
    )
    replacements = {
        "ICHARG": "11",
        "NCORE": "1",
        "LWANNIER90": "True",
        "LWRITE_MMN_AMN": "True",
        "NUM_WANN": str(num_wann),
        "NBANDS": str(nbands),
    }
    seen = set()
    messages = []
    original = incar_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    lines = []
    index = 0
    while index < len(original):
        raw_line = original[index]
        stripped = raw_line.strip()
        active = stripped and not stripped.startswith(("#", "!")) and "=" in raw_line
        if active:
            key = raw_line.split("=", 1)[0].strip().upper()
            if key == "WANNIER90_WIN":
                index = strip_wannier90_win(original, index)
                continue
            if key == "NPAR":
                messages.append("removed active NPAR in INCAR for projection stage")
                index += 1
                continue
            if key == "ISPIN" and spinor:
                messages.append("removed active ISPIN because LSORBIT/LNONCOLLINEAR is enabled")
                index += 1
                continue
            if key in replacements:
                lines.append(f"{key:<20} = {replacements[key]}")
                seen.add(key)
                index += 1
                continue
        lines.append(raw_line)
        index += 1

    lines.append("")
    lines.append("# Wannier90 projection settings added by auto-MLWFs skill")
    for key, value in replacements.items():
        if key not in seen:
            lines.append(f"{key:<20} = {value}")
    lines.extend(make_wannier_win(projections, kmesh_tol))
    incar_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return messages


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Calculation root containing scf/ and proj/.")
    parser.add_argument("--scf-dir", default="scf")
    parser.add_argument("--proj-dir", default="proj")
    parser.add_argument(
        "--run-np",
        type=int,
        required=True,
        help="MPI ranks that will actually run VASP; used to round NBANDS.",
    )
    parser.add_argument("--projections", default=None, help="Projection list; defaults to scf/outputs/projections.txt.")
    parser.add_argument(
        "--existing",
        choices=["agent", "fail", "backup", "overwrite", "continue"],
        default="agent",
        help="On collision, default to <proj-dir>_agent; fail if both paths exist.",
    )
    parser.add_argument("--kmesh-tol", type=float, default=0.0001)
    parser.add_argument("--dry-run", action="store_true", help="Compute settings without copying or editing files.")
    args = parser.parse_args()
    if args.run_np <= 0:
        raise ValueError("--run-np must be positive")

    root = Path(args.root).resolve()
    scf_dir = root / args.scf_dir
    requested_proj_dir = root / args.proj_dir
    proj_dir = agent_output_dir(requested_proj_dir) if args.existing == "agent" else requested_proj_dir
    existing_mode = "fail" if args.existing == "agent" else args.existing
    projection_path = Path(args.projections).resolve() if args.projections else scf_dir / "outputs" / "projections.txt"

    projections = read_projections(projection_path)
    atom_counts = parse_poscar_counts(scf_dir / "POSCAR")
    flags = incar_flags(scf_dir / "INCAR")
    spinor = flags["LSORBIT"] or flags["LNONCOLLINEAR"]
    num_wann = compute_num_wann(projections, atom_counts, soc=spinor)
    nelect = read_nelect(scf_dir / "OUTCAR")
    nbands = nbands_from_nelect(nelect, args.run_np, lower_bound=num_wann)

    if args.dry_run:
        print(f"proj_dir: {proj_dir}")
        print(f"projections: {', '.join(projections)}")
        print(f"nelect: {nelect:g}")
        print(f"num_wann: {num_wann}")
        print(f"spinor: {str(spinor).lower()}")
        print(f"nbands: {nbands}")
        print(f"run_np: {args.run_np}")
        print(f"vasp_mode: {'ncl' if spinor else 'std'}")
        print(f"kmesh_tol: {args.kmesh_tol:g}")
        print("dry_run: true")
        return

    copy_scf_to_proj(scf_dir, proj_dir, existing_mode)
    messages = update_incar(proj_dir / "INCAR", projections, num_wann, nbands, args.kmesh_tol)
    mode = vasp_mode(proj_dir / "INCAR")
    summary_path = write_stage_summary(
        proj_dir,
        {
            "schema_version": 1,
            "status": "prepared",
            "stage": "projection",
            "requested_workdir": str(requested_proj_dir),
            "workdir": str(proj_dir),
            "used_agent_suffix": proj_dir != requested_proj_dir,
            "projections": projections,
            "nelect": nelect,
            "num_wann": num_wann,
            "spinor": spinor,
            "nbands": nbands,
            "run_np": args.run_np,
            "vasp_mode": mode,
            "kmesh_tol": args.kmesh_tol,
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
            "notes": messages,
        },
    )
    print(f"proj_dir: {proj_dir}")
    print(f"projections: {', '.join(projections)}")
    print(f"nelect: {nelect:g}")
    print(f"num_wann: {num_wann}")
    print(f"spinor: {str(spinor).lower()}")
    print(f"nbands: {nbands}")
    print(f"run_np: {args.run_np}")
    print(f"vasp_mode: {mode}")
    print(f"kmesh_tol: {args.kmesh_tol:g}")
    for message in messages:
        print(f"note: {message}")
    print("scheduler_script: user_or_agent_managed")
    print(f"stage_summary: {summary_path}")


if __name__ == "__main__":
    main()
