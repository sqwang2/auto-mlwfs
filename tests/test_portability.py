"""Exercise an installed skill outside both the repository and calculation."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml


SKILL = Path(os.environ.get(
    "AUTO_MLWFS_TEST_SKILL_DIR",
    Path(__file__).resolve().parents[1] / "skills" / "auto-mlwfs",
)).resolve()
sys.path.insert(0, str(SKILL / "scripts"))

from dos_and_analysis import load_config
from optimize_windows import load_yaml, optimization_paths
from workflow_io import resolve_config


class InstalledSkillTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="auto mlwfs test ")
        self.addCleanup(self.temporary.cleanup)
        base = Path(self.temporary.name)
        self.installed = base / "installed skill"
        shutil.copytree(SKILL, self.installed, ignore=shutil.ignore_patterns("__pycache__"))
        self.root = base / "calculation"
        self.cwd = base / "unrelated directory"
        self.cwd.mkdir()
        (self.root / "band").mkdir(parents=True)
        (self.root / "wann").mkdir()

    def run_script(self, name, *args):
        return subprocess.run(
            [sys.executable, str(self.installed / "scripts" / name), *args],
            cwd=self.cwd,
            env={**os.environ, "MPLBACKEND": "Agg"},
            text=True, capture_output=True, timeout=90,
        )

    def write_bands(self, scale=1.0):
        energies = [-4.5, -3.2, -2.1, -1.3, -0.6, -0.1, 0.4, 1.1, 1.9, 2.8, 4.0, 5.3]
        eigenval = ["synthetic test data"] * 5 + ["12 2 12"]
        for k in range(2):
            eigenval.extend(["", f"{k} 0 0 1"])
            eigenval.extend(f"{i} {e + k * 0.03} 0" for i, e in enumerate(energies, 1))
        (self.root / "band" / "EIGENVAL").write_text("\n".join(eigenval) + "\n")
        (self.root / "band" / "OUTCAR").write_text(" E-fermi : 0.0000\n")
        blocks = [f"0 {e * scale}\n1 {(e + 0.03) * scale}" for e in energies]
        (self.root / "wann" / "wannier90_band.dat").write_text("\n\n".join(blocks) + "\n")
        (self.root / "wann" / "window_summary.json").write_text(json.dumps({
            "dis_froz_min": -2.0, "dis_froz_max": 2.0,
        }))

    def test_all_cli_entrypoints_load_after_relocation(self):
        for name in [
            "check_success.py", "dos_and_analysis.py", "lesa_fitting.py",
            "optimize_windows.py", "plot_fitting.py", "prepare_projection.py",
            "prepare_wannier.py", "prepare_band_from_wannier_kpt.py", "vasp_band.py",
        ]:
            with self.subTest(script=name):
                result = self.run_script(name, "--help")
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_lesa_uses_installed_defaults_from_unrelated_cwd(self):
        self.write_bands()
        result = self.run_script("lesa_fitting.py", "--root", str(self.root))
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads((self.root / "outputs" / "fitting_lesa_summary.json").read_text())
        self.assertEqual(summary["result"]["status"], "pass")
        self.assertEqual(summary["frozen_window"]["ratio"], 1.0)
        self.assertFalse((self.cwd / "outputs").exists())

    def test_lesa_rejects_mismatched_bands(self):
        self.write_bands(scale=1.5)
        result = self.run_script("lesa_fitting.py", "--root", str(self.root))
        self.assertEqual(result.returncode, 1, result.stderr)
        summary = json.loads((self.root / "outputs" / "fitting_lesa_summary.json").read_text())
        self.assertEqual(summary["result"]["status"], "fail")

    def test_custom_lesa_config_is_relative_to_calculation_root(self):
        self.write_bands()
        config = yaml.safe_load((self.installed / "configuration" / "fitting.yaml").read_text())
        config["paths"]["output_dir"] = "custom results"
        (self.root / "custom.yaml").write_text(yaml.safe_dump(config))
        result = self.run_script("lesa_fitting.py", "--root", str(self.root), "--config", "custom.yaml")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.root / "custom results" / "fitting_lesa_summary.json").exists())
        self.assertFalse((self.root / "outputs").exists())

    def test_optimization_fitting_defaults_and_overrides(self):
        config = load_yaml(SKILL / "configuration" / "optimization.yaml")
        defaults = optimization_paths(self.root, config)
        self.assertTrue(defaults["fitting_config"].is_file())
        self.assertEqual(defaults["band"], self.root / "band")
        config["paths"]["fitting_config"] = "custom.yaml"
        self.assertEqual(optimization_paths(self.root, config)["fitting_config"], self.root / "custom.yaml")
        absolute = self.cwd / "external.yaml"
        config["paths"]["fitting_config"] = str(absolute)
        self.assertEqual(optimization_paths(self.root, config)["fitting_config"], absolute)

    def test_projection_config_missing_override_fails(self):
        bundled = resolve_config(self.root, None, "projection.yaml")
        self.assertEqual(load_config(bundled)["analysis"]["mode"], "dos_manifold")
        with self.assertRaises(FileNotFoundError):
            load_config(resolve_config(self.root, "missing.yaml", "projection.yaml"))


if __name__ == "__main__":
    unittest.main()
