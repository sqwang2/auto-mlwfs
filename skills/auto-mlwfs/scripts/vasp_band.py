#!/usr/bin/env python3
"""Prepare a VASP band-structure calculation from an SCF directory."""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
from pathlib import Path

from workflow_io import agent_output_dir, write_stage_summary


VASP_OUTPUT_PATTERNS = [
    "WAVECAR",
    "vasprun.xml",
    "OUTCAR",
    "OSZICAR",
    "DOSCAR",
    "PROCAR",
    "EIGENVAL",
    "REPORT",
    "PCDAT",
    "XDATCAR",
    "IBZKPT",
    "slurm-*.out",
    "*.tmp",
]

KPATH_TASKS = {
    "1d": "301",
    "2d": "302",
    "bulk": "303",
    "3d": "303",
}


def read_nelect(outcar_path: Path) -> float:
    if not outcar_path.exists():
        raise FileNotFoundError(f"OUTCAR not found for NELECT detection: {outcar_path}")
    text = outcar_path.read_text(encoding="utf-8", errors="ignore")
    matches = re.findall(r"\bNELECT\s*=\s*([0-9]+(?:\.[0-9]*)?)", text)
    if not matches:
        raise ValueError(f"Could not find NELECT in {outcar_path}")
    return float(matches[-1])


def round_up_to_multiple(value: float, multiple: int) -> int:
    if multiple <= 0:
        return int(math.ceil(value))
    return int(math.ceil(value / multiple) * multiple)


def nbands_from_nelect(nelect: float, run_np: int, factor: float = 2.0, lower_bound: int = 0) -> int:
    target = max(float(lower_bound), factor * nelect)
    return round_up_to_multiple(target, int(run_np))


def incar_flags(incar_path: Path) -> dict:
    flags = {"LSORBIT": False, "LNONCOLLINEAR": False}
    if not incar_path.exists():
        return flags
    for raw_line in incar_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.split("#", 1)[0].split("!", 1)[0]
        if "=" not in line:
            continue
        key, value = [part.strip().upper() for part in line.split("=", 1)]
        if key in flags:
            flags[key] = value in {"T", ".TRUE.", "TRUE", "1", "YES"}
    return flags


def vasp_mode(incar_path: Path) -> str:
    flags = incar_flags(incar_path)
    if flags["LSORBIT"] or flags["LNONCOLLINEAR"]:
        return "ncl"
    return "std"


def copy_scf_to_band(scf_dir: Path, band_dir: Path, mode: str) -> None:
    if not scf_dir.exists():
        raise FileNotFoundError(f"SCF directory not found: {scf_dir}")
    if band_dir.exists():
        if mode == "fail":
            raise FileExistsError(f"Band directory already exists: {band_dir}")
        if mode == "backup":
            index = 1
            while True:
                backup = band_dir.with_name(f"{band_dir.name}.bak{index}")
                if not backup.exists():
                    band_dir.rename(backup)
                    break
                index += 1
        elif mode == "overwrite":
            shutil.rmtree(band_dir)
        elif mode == "continue":
            return
    shutil.copytree(scf_dir, band_dir)


def clean_band_dir(band_dir: Path) -> None:
    for pattern in VASP_OUTPUT_PATTERNS:
        for path in band_dir.glob(pattern):
            if path.is_file() or path.is_symlink():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)


def update_incar(incar_path: Path, nbands: int | None = None) -> None:
    if not incar_path.exists():
        raise FileNotFoundError(f"INCAR not found: {incar_path}")
    spinor = any(incar_flags(incar_path).values())
    seen_icharg = False
    seen_nedos = False
    seen_nbands = False
    lines = []
    for raw_line in incar_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = raw_line.strip()
        active = stripped and not stripped.startswith(("#", "!")) and "=" in raw_line
        if active:
            key = raw_line.split("=", 1)[0].strip().upper()
            if key == "ICHARG":
                lines.append(f"{key:<20} = 11")
                seen_icharg = True
                continue
            if key == "NBANDS" and nbands is not None:
                lines.append(f"{key:<20} = {nbands}")
                seen_nbands = True
                continue
            if key == "ISPIN" and spinor:
                continue
            if key == "NEDOS":
                seen_nedos = True
        lines.append(raw_line)
    if (
        not seen_icharg
        or not seen_nedos
        or (nbands is not None and not seen_nbands)
    ):
        lines.append("")
        lines.append("# Band calculation settings added by auto-MLWFs skill")
        if not seen_icharg:
            lines.append(f"{'ICHARG':<20} = 11")
        if nbands is not None and not seen_nbands:
            lines.append(f"{'NBANDS':<20} = {nbands}")
        if not seen_nedos:
            lines.append(f"{'NEDOS':<20} = 3001")
    incar_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_vaspkit_kpath(band_dir: Path, task: str, vaspkit_command: str) -> None:
    completed = subprocess.run(
        [vaspkit_command],
        input=f"{task}\n",
        cwd=band_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (band_dir / "vaspkit-kpath.log").write_text(completed.stdout, encoding="utf-8", errors="ignore")
    if completed.returncode != 0:
        raise RuntimeError(f"VASPKIT failed with exit code {completed.returncode}; see vaspkit-kpath.log")
    kpath = band_dir / "KPATH.in"
    if not kpath.exists():
        raise FileNotFoundError("VASPKIT did not generate KPATH.in")
    shutil.copy2(kpath, band_dir / "KPOINTS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".", help="Calculation root containing scf/ and band/.")
    parser.add_argument("--scf-dir", default="scf")
    parser.add_argument("--band-dir", default="band")
    parser.add_argument(
        "--run-np",
        type=int,
        required=True,
        help="MPI ranks that will actually run VASP; used to round NBANDS.",
    )
    parser.add_argument("--dimension", choices=sorted(KPATH_TASKS), default="bulk")
    parser.add_argument("--vaspkit-task", default=None, help="Override VASPKIT task number, e.g. 303.")
    parser.add_argument("--vaspkit-command", default="vaspkit")
    parser.add_argument(
        "--existing",
        choices=["agent", "fail", "backup", "overwrite", "continue"],
        default="agent",
        help="On collision, default to <band-dir>_agent; fail if both paths exist.",
    )
    parser.add_argument("--skip-vaspkit", action="store_true", help="Prepare inputs without running VASPKIT.")
    args = parser.parse_args()
    if args.run_np <= 0:
        raise ValueError("--run-np must be positive")

    root = Path(args.root).resolve()
    scf_dir = root / args.scf_dir
    requested_band_dir = root / args.band_dir
    band_dir = agent_output_dir(requested_band_dir) if args.existing == "agent" else requested_band_dir
    existing_mode = "fail" if args.existing == "agent" else args.existing
    task = args.vaspkit_task or KPATH_TASKS[args.dimension]
    nelect = read_nelect(scf_dir / "OUTCAR")
    nbands = nbands_from_nelect(nelect, args.run_np)

    copy_scf_to_band(scf_dir, band_dir, existing_mode)
    clean_band_dir(band_dir)
    update_incar(band_dir / "INCAR", nbands=nbands)
    if not args.skip_vaspkit:
        run_vaspkit_kpath(band_dir, task, args.vaspkit_command)
    summary_path = write_stage_summary(
        band_dir,
        {
            "schema_version": 1,
            "status": "prepared",
            "stage": "vasp-band-path",
            "requested_workdir": str(requested_band_dir),
            "workdir": str(band_dir),
            "used_agent_suffix": band_dir != requested_band_dir,
            "kpath_task": task,
            "nelect": nelect,
            "nbands": nbands,
            "run_np": args.run_np,
            "vasp_mode": vasp_mode(band_dir / "INCAR"),
            "scheduler_required": False,
            "next_action": "prepare_projection",
        },
    )

    print(f"band_dir: {band_dir}")
    print(f"kpath_task: {task}")
    print(f"nelect: {nelect:g}")
    print(f"nbands: {nbands}")
    print(f"run_np: {args.run_np}")
    print(f"vasp_mode: {vasp_mode(band_dir / 'INCAR')}")
    print("scheduler_script: not_required_for_path_preparation")
    print(f"stage_summary: {summary_path}")


if __name__ == "__main__":
    main()
