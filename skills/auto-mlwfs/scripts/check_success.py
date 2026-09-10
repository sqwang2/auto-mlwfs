#!/usr/bin/env python3
"""Check whether a workflow stage finished successfully."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


VASP_BAND_TIMING_FOOTER = "General timing and accounting informations for this job"
WANNIER_KPOINT_NEIGHBOURS = "Finished setting up k-point neighbours."
WANNIER_SETUP_EXIT = "Exiting wannier_setup in wannier90"
WANNIER_ALL_DONE = "All done: wannier90 exiting"


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def active_incar_value(incar_path: Path, key: str) -> str | None:
    target = key.upper()
    for raw_line in read_text(incar_path).splitlines():
        line = raw_line.split("#", 1)[0].split("!", 1)[0].strip()
        if not line or "=" not in line:
            continue
        lhs, rhs = line.split("=", 1)
        if lhs.strip().upper() == target:
            return rhs.strip().split()[0] if rhs.strip() else ""
    return None


def check_vasp_band(workdir: Path) -> tuple[bool, list[str]]:
    messages: list[str] = []

    outcar = workdir / "OUTCAR"
    incar = workdir / "INCAR"
    kpoints = workdir / "KPOINTS"

    if not outcar.exists():
        messages.append("FAIL OUTCAR is missing.")
    elif VASP_BAND_TIMING_FOOTER not in read_text(outcar):
        messages.append(f"FAIL OUTCAR does not contain '{VASP_BAND_TIMING_FOOTER}'.")
    else:
        messages.append("PASS OUTCAR contains VASP timing footer.")

    if not incar.exists():
        messages.append("FAIL INCAR is missing.")
    else:
        icharg = active_incar_value(incar, "ICHARG")
        if icharg == "11":
            messages.append("PASS INCAR has active ICHARG = 11.")
        else:
            messages.append(f"FAIL INCAR active ICHARG is {icharg!r}, expected '11'.")

    if not kpoints.exists():
        messages.append("FAIL KPOINTS is missing.")
    else:
        kpoints_text = read_text(kpoints)
        if re.search(r"^\s*Line-Mode\s*$", kpoints_text, re.IGNORECASE | re.MULTILINE):
            messages.append("PASS KPOINTS is Line-Mode.")
        elif re.search(r"^\s*(Reciprocal|Cartesian)\s*$", kpoints_text, re.IGNORECASE | re.MULTILINE):
            messages.append("PASS KPOINTS is an explicit coordinate list.")
        else:
            messages.append("FAIL KPOINTS is neither Line-Mode nor an explicit coordinate list.")

    ok = all(message.startswith("PASS") for message in messages)
    return ok, messages


def check_projection(workdir: Path, seedname: str = "wannier90") -> tuple[bool, list[str]]:
    messages: list[str] = []

    outcar = workdir / "OUTCAR"
    wout = workdir / f"{seedname}.wout"

    if not outcar.exists():
        messages.append("FAIL OUTCAR is missing.")
    elif VASP_BAND_TIMING_FOOTER not in read_text(outcar):
        messages.append(f"FAIL OUTCAR does not contain '{VASP_BAND_TIMING_FOOTER}'.")
    else:
        messages.append("PASS OUTCAR contains VASP timing footer.")

    if not wout.exists():
        messages.append(f"FAIL {seedname}.wout is missing.")
    else:
        wout_text = read_text(wout)
        if WANNIER_KPOINT_NEIGHBOURS in wout_text:
            messages.append(f"PASS {seedname}.wout contains Wannier k-point neighbour setup completion.")
        else:
            messages.append(f"FAIL {seedname}.wout does not contain '{WANNIER_KPOINT_NEIGHBOURS}'.")
        if WANNIER_SETUP_EXIT in wout_text:
            messages.append(f"PASS {seedname}.wout contains Wannier setup exit marker.")
        else:
            messages.append(f"FAIL {seedname}.wout does not contain '{WANNIER_SETUP_EXIT}'.")

    for suffix in ("amn", "mmn", "eig"):
        path = workdir / f"{seedname}.{suffix}"
        if path.exists() and path.stat().st_size > 0:
            messages.append(f"PASS {seedname}.{suffix} exists and is non-empty.")
        elif path.exists():
            messages.append(f"FAIL {seedname}.{suffix} exists but is empty.")
        else:
            messages.append(f"FAIL {seedname}.{suffix} is missing.")

    ok = all(message.startswith("PASS") for message in messages)
    return ok, messages


def check_wannier(workdir: Path, seedname: str = "wannier90") -> tuple[bool, list[str]]:
    messages: list[str] = []
    wout = workdir / f"{seedname}.wout"
    band_file = workdir / f"{seedname}_band.dat"

    if not wout.exists():
        messages.append(f"FAIL {seedname}.wout is missing.")
    elif WANNIER_ALL_DONE not in read_text(wout):
        messages.append(f"FAIL {seedname}.wout does not contain '{WANNIER_ALL_DONE}'.")
    else:
        messages.append(f"PASS {seedname}.wout contains the Wannier90 completion marker.")

    if band_file.exists() and band_file.stat().st_size > 0:
        messages.append(f"PASS {band_file.name} exists and is non-empty.")
    elif band_file.exists():
        messages.append(f"FAIL {band_file.name} exists but is empty.")
    else:
        messages.append(f"FAIL {band_file.name} is missing.")

    ok = all(message.startswith("PASS") for message in messages)
    return ok, messages


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["vasp-band", "band", "projection", "wannier"], required=True)
    parser.add_argument("--workdir", default=".", help="Stage calculation directory, e.g. band/")
    parser.add_argument("--seedname", default="wannier90", help="Wannier90 seedname for projection checks.")
    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    if not workdir.exists():
        print(f"ERROR workdir does not exist: {workdir}")
        return 2

    if args.stage in {"vasp-band", "band"}:
        stage_name = "vasp-band"
        ok, messages = check_vasp_band(workdir)
    elif args.stage == "projection":
        stage_name = "projection"
        ok, messages = check_projection(workdir, args.seedname)
    else:
        stage_name = "wannier"
        ok, messages = check_wannier(workdir, args.seedname)

    print(("PASS" if ok else "FAIL") + f" {stage_name}")
    for message in messages:
        print(f"- {message}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
