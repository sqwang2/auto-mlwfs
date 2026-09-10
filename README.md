# Auto MLWFs

An agent skill for VASP-to-Wannier90 workflows: PDOS-based projection selection,
input preparation, band comparison, LESA fitting assessment, and Bayesian
optimization of disentanglement and frozen windows.

Auto MLWFs 是面向 Codex、Claude Code 等智能体的科研技能，辅助完成
VASP 到 Wannier90 的工作流、拟合评估与能窗优化。

## Install / 一行安装

With Node.js/npm and Git installed, run:

```bash
npx skills add sqwang2/auto-mlwfs --skill auto-mlwfs -a codex -g -y
```

For Claude Code:

```bash
npx skills add sqwang2/auto-mlwfs --skill auto-mlwfs -a claude-code -g -y
```

Omit `-g` to install into the current project. To choose an agent interactively:

```bash
npx skills add sqwang2/auto-mlwfs --skill auto-mlwfs
```

The third-party [skills CLI](https://github.com/vercel-labs/skills) installs the
skill instructions, scripts, and configuration. It does not install scientific
software or Python dependencies. Start a new agent session if the skill does
not appear immediately.

这条命令安装技能本身。完整计算还需要下方的 Python 环境和外部科研软件。

## Runtime / 计算环境

- Python 3.11 or newer; Python dependencies are listed in the skill's
  [requirements.txt](skills/auto-mlwfs/requirements.txt).
- An existing VASP SCF calculation with the outputs needed by each stage.
- A licensed VASP installation, Wannier90, and VASPKIT for automated k-path
  generation. VASPKIT can be skipped when a valid path is supplied manually.
- MPI and Slurm/PBS as required by the execution environment.

For a default global Codex installation on macOS/Linux, create an isolated
Python environment on the machine that will run the scripts:

```bash
export AUTO_MLWFS_SKILL_DIR="$HOME/.agents/skills/auto-mlwfs"
python3 -m venv "$HOME/.venvs/auto-mlwfs"
source "$HOME/.venvs/auto-mlwfs/bin/activate"
python -m pip install -r "$AUTO_MLWFS_SKILL_DIR/requirements.txt"
```

The skills CLI uses the shared `~/.agents/skills/` directory for Codex global
installs. For another agent, an older installer, or a custom installation, set
`AUTO_MLWFS_SKILL_DIR` to the actual directory reported by the installer that
contains this skill's `SKILL.md`. Use a Python 3.11+ command
in place of `python3` if your default Python is older. These shell examples
target macOS/Linux and Linux HPC systems.

安装目录与计算目录相互独立。集群运行时，需要在执行节点可访问的位置准备技能
和 Python 环境，并在作业脚本中使用实际的解释器与脚本路径。

VASP binaries, POTCAR files, other external executables, and cluster credentials
are not distributed by this project. Users supply software under its applicable
license and configure their own executable paths, queues, and resources.

## Use / 使用

Ask your agent, for example:

> 使用 auto-mlwfs，分析当前项目 scf/ 的 PDOS 并推荐 Wannier 投影轨道。
> 先检查环境与输入文件，再准备后续计算。

Or run a script directly:

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/dos_and_analysis.py" --workdir /path/to/project/scf
```

The [skill instructions](skills/auto-mlwfs/SKILL.md) describe the full staged
workflow and [scheduler setup](skills/auto-mlwfs/docs/scheduler_guidance.md).
The agent prepares jobs using your environment and confirms submission with you.
Scientific outputs require inspection; successful installation is not a
validation of a material's Wannier representation.

## Configuration

Bundled defaults are resolved relative to the installed skill, independently of
the current working directory. Copy a YAML file into your calculation project
to customize it, and supply `--config` explicitly.

- LESA and optimization resolve relative `--config` paths against `--root`.
- DOS resolves relative `--config` paths against the shell's current directory.
- Calculation paths inside YAML are relative to `--root` or DOS `--workdir`.
- Optimization `paths.fitting_config: null` uses the bundled LESA configuration;
  a custom path is relative to `--root` unless absolute.
- Use the same custom optimization configuration for `propose`, `run-trial`,
  `collect`, and `status`.

The [AgI phase-2 example](skills/auto-mlwfs/examples/agi/optimization_phase2.yaml)
references a historical trial. Replace its baseline path before using it on
your own data. Default Bayesian optimization has no evaluation limit; set
`bo.max_evaluations` in your project configuration when a finite budget is needed.

## Development

```bash
git clone https://github.com/sqwang2/auto-mlwfs.git
cd auto-mlwfs
python -m pip install -r skills/auto-mlwfs/requirements.txt
python -m unittest discover -s tests -v
npx skills add . --list
```

The tests exercise installation portability and synthetic analysis data without
submitting VASP or Wannier90 jobs. They do not replace scientific validation
on real materials or testing on a particular cluster.

## License and Citation

[MIT](LICENSE). External scientific software retains its own license.
Use [CITATION.cff](CITATION.cff) to cite this software; cite the underlying
scientific methods and external software separately where applicable.
