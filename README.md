# Auto-MLWFs Skill

**简体中文** | [English](README.en.md)

[![Validate Skill](https://github.com/sqwang2/auto-mlwfs/actions/workflows/ci.yml/badge.svg)](https://github.com/sqwang2/auto-mlwfs/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/Python-3.11%2B-3776AB?logo=python&logoColor=white)](skills/auto-mlwfs/requirements.txt)
[![License: MIT](https://img.shields.io/badge/License-MIT-2A7F62)](LICENSE)

## 安装与接入

准备好 Node.js/npm 和 Git 后，一行安装到 **Codex**：

```bash
npx skills add sqwang2/auto-mlwfs --skill auto-mlwfs -a codex -g -y
```

安装到 **Claude Code**：

```bash
npx skills add sqwang2/auto-mlwfs --skill auto-mlwfs -a claude-code -g -y
```

`-g` 表示用户级安装；去掉它即可安装到当前项目。安装后，在 agent 中请求使用 `auto-mlwfs`。若技能没有立即显示，可开启新会话。

<details>
<summary><strong>卡在 Cloning repository？使用 GitHub 官方源码包安装</strong></summary>

如果终端无法连接 `github.com`，但可以访问 `codeload.github.com`，先按 `Ctrl+C` 结束卡住的命令，再执行：

```bash
npx --yes skills@latest add https://codeload.github.com/sqwang2/auto-mlwfs/tar.gz/refs/heads/main --skill auto-mlwfs -a codex -g -y
```

这会从 GitHub 官方下载源码包并安装，绕过 Git 克隆。它不是第三方镜像；更换 npm 镜像源不会解决 GitHub 克隆阶段的连接问题。安装到 Claude Code 时，将 `-a codex` 改为 `-a claude-code`。

</details>

安装由第三方 [skills CLI](https://github.com/vercel-labs/skills) 完成，包含技能指令、脚本和配置。Python 依赖及 VASP、Wannier90 等外部程序需要按下方[环境要求](#环境要求)单独准备。

## 项目简介

面向 AI agent 的 VASP–Wannier90 最大局域化 Wannier 函数（MLWF）工作流。

本项目把投影轨道选择、VASP/Wannier90 输入准备、能窗构造、拟合评估和贝叶斯优化组织成一套可复用的 Skill。它适合从已经完成的 VASP SCF 计算出发，在 agent 与用户共同确认计算资源和提交权限的前提下，完成后续 Wannier 化流程。

> 当前状态：研究与工程验证阶段。建议先在已有算例上测试，再用于新的材料体系。

## 主要功能

- 根据 `vasprun.xml` 的 DOS/PDOS 自动推荐 Wannier 投影轨道；
- 支持普通标量、自旋轨道耦合（SOC）和非共线计算的输入准备；
- 自动准备 VASP band、projection 和 Wannier90 输入；
- 根据 `NELECT`、`NUM_WANN` 和实际 MPI ranks 设置 `NBANDS`；
- 从投影能带构造物理初始 disentanglement/frozen window；
- 生成 DFT–Wannier 拟合图，并通过 LESA 给出机器可读的拟合判断；
- 在初始结果不合格时，对四个能窗边界进行贝叶斯优化（BO）；
- 支持 agent 根据当前环境适配 Slurm 或 PBS，而不是依赖固定的提交模板；
- 通过 `stage_summary.json` 和 BO batch JSON 向上层 agent 返回结构化结果。

## 工作流概览

```text
completed VASP SCF
        |
        +--> DOS/PDOS analysis --> projection recommendation
        |
        +--> VASP line-mode k-path preparation
        |
        +--> VASP projection calculation
        |
        +--> Wannier90 input and physics-window preparation
        |
        +--> Wannier90 calculation
        |
        +--> final VASP band calculation on Wannier90 k-points
        |
        +--> fitting plot + LESA assessment
                         |
                         +--> pass: deliver results
                         |
                         `--> fail: four-dimensional BO
```

## 目录结构

```text
auto-mlwfs/
├── README.md                        # 中文说明
├── README.en.md                     # English documentation
├── LICENSE
├── CITATION.cff
├── skills/
│   └── auto-mlwfs/
│       ├── SKILL.md                 # agent 的核心工作指令
│       ├── LICENSE
│       ├── requirements.txt
│       ├── configuration/
│       │   ├── projection.yaml      # 投影轨道选择参数
│       │   ├── fitting.yaml         # LESA 拟合判据
│       │   └── optimization.yaml    # 默认 BO 配置
│       ├── docs/                    # 工作流、判据与故障处理
│       ├── examples/agi/
│       │   └── optimization_phase2.yaml
│       └── scripts/
│           ├── dos_and_analysis.py
│           ├── vasp_band.py
│           ├── prepare_projection.py
│           ├── prepare_wannier.py
│           ├── prepare_band_from_wannier_kpt.py
│           ├── plot_fitting.py
│           ├── lesa_fitting.py
│           ├── optimize_windows.py
│           ├── check_success.py
│           └── ...                 # 共享解析与输出辅助模块
└── tests/                           # 安装可移植性与合成数据测试
```

## 环境要求

### 外部程序

- 具有相应使用许可的 VASP；
- VASPKIT，用于自动生成高对称 k 路径；已有有效路径时可跳过；
- Wannier90；
- 可选的 Slurm 或 PBS 调度环境。

程序位置、启动命令、队列、核数、内存和环境模块均由用户在实际运行环境中指定，本项目不内置站点专用路径。

VASP 可执行文件、POTCAR、其他外部程序和集群凭据不随本项目分发。

### Python 依赖

使用 **Python 3.11 或更新版本**。基础分析与绘图需要：

```text
numpy >= 2.0
pymatgen
matplotlib
PyYAML
```

运行贝叶斯优化还需要：

```text
scikit-optimize
```

下面以 macOS/Linux 上通过 skills CLI 用户级安装到 Codex 为例，创建独立的 Python 环境：

```bash
export AUTO_MLWFS_SKILL_DIR="$HOME/.agents/skills/auto-mlwfs"
python3 -m venv "$HOME/.venvs/auto-mlwfs"
source "$HOME/.venvs/auto-mlwfs/bin/activate"
python -m pip install -r "$AUTO_MLWFS_SKILL_DIR/requirements.txt"
```

请确认 `python3` 对应 Python 3.11+。其他 agent、项目级安装或旧版安装工具的目录可能不同，`AUTO_MLWFS_SKILL_DIR` 应指向安装工具实际报告的、包含 `SKILL.md` 的目录。

完整依赖版本范围见 [`requirements.txt`](skills/auto-mlwfs/requirements.txt)。在集群执行时，技能目录和 Python 环境需要对执行节点可见；作业脚本应使用对应解释器。

## 输入要求

典型项目从一个已完成的 `scf/` 目录开始。根据所执行的阶段，通常需要：

```text
scf/
├── INCAR
├── POSCAR
├── KPOINTS
├── POTCAR                           # 用户按许可自行提供
├── OUTCAR
├── vasprun.xml
├── EIGENVAL
└── CHGCAR
```

DOS/PDOS 分析要求 `vasprun.xml` 中包含实际的投影态密度数据。

技能安装目录与计算目录相互独立。下面的命令从计算项目根目录执行，使用环境设置中定义的 `AUTO_MLWFS_SKILL_DIR`：

```text
project/
└── scf/
```

例如，可以向 agent 发出请求：“使用 auto-mlwfs，分析当前项目 scf/ 的 PDOS 并推荐投影轨道；检查输入与环境后，再准备后续计算。”

## 快速开始

### 1. 分析 DOS/PDOS 并推荐投影轨道

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/dos_and_analysis.py" \
  --workdir scf \
  --config "$AUTO_MLWFS_SKILL_DIR/configuration/projection.yaml"
```

主要输出位于 `scf/outputs/`：

- `projection_summary.json`
- `projections.txt`
- `dos_projected.csv`
- `dos_projected_weights.csv`
- `DOS.png`

### 2. 准备高对称 k 路径

先询问用户 VASP 实际使用的 MPI ranks：

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/vasp_band.py" --root . --run-np <actual-mpi-ranks>
```

该阶段使用 VASPKIT 生成 line-mode k 路径，供后续 Wannier90 使用；它不是最终用于拟合的 VASP band 计算。

### 3. 准备并运行 projection 计算

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_projection.py" \
  --root . \
  --run-np <actual-mpi-ranks>
```

脚本读取 `scf/outputs/projections.txt`，计算 `NUM_WANN`，并准备 `proj/`。随后由 agent 检测 Slurm/PBS、询问缺失的执行参数、生成提交脚本，并在用户明确确认后提交。

完成后检查：

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/check_success.py" --stage projection --workdir proj
```

### 4. 准备并运行 Wannier90

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_wannier.py" --root .
```

该脚本会根据投影能带生成物理初始能窗，写入 `wann/wannier90.win`，并生成 `wann/window_summary.json`。随后按当前集群环境运行 Wannier90。

完成后检查：

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/check_success.py" --stage wannier --workdir wann
```

### 5. 准备最终 VASP band 计算

Wannier90 生成 `wannier90_band.kpt` 后执行：

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/prepare_band_from_wannier_kpt.py" \
  --root . \
  --wann-dir wann \
  --band-dir band \
  --run-np <actual-mpi-ranks>
```

脚本以 Wannier90 实际采样的 k 点重写 `band/KPOINTS`。之后再由 agent 与用户确认资源和提交权限，运行最终 VASP band 任务。

完成后检查：

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/check_success.py" --stage vasp-band --workdir band
```

### 6. 绘图并进行 LESA 评估

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/plot_fitting.py" --root . --vasp-dir band --wann-dir wann

python "$AUTO_MLWFS_SKILL_DIR/scripts/lesa_fitting.py" \
  --root . \
  --config "$AUTO_MLWFS_SKILL_DIR/configuration/fitting.yaml"
```

主要输出包括：

- `<project-name>_fitting.png`
- `outputs/fitting_lesa_summary.json`
- `outputs/fitting_lesa_report.txt`

拟合图用于人工检查；LESA JSON 是 agent 的默认机器可读判断。两者应配套保留。

当前拟合解析不能完整评估共线自旋极化 `ISPIN=2` 的两个自旋通道，因此不应把输入准备支持范围视为所有自旋模式的端到端验证承诺。

## NBANDS 规则

`NBANDS` 不从旧 INCAR 直接沿用，而是从 `scf/OUTCAR` 读取 `NELECT` 后重新计算。

band 和最终 band 阶段：

```text
NBANDS = ceil_to_multiple(2 × NELECT, actual_mpi_ranks)
```

projection 阶段：

```text
NBANDS = ceil_to_multiple(max(2 × NELECT, NUM_WANN), actual_mpi_ranks)
```

这里使用的是程序实际运行的 MPI ranks，而不是调度器申请的总核数。如果因内存问题改变实际 ranks，应重新执行相应的准备脚本。

## 调度器与权限边界

本项目不会强行绑定某台集群的 submit YAML。Agent 应当：

1. 检测当前环境使用 Slurm 还是 PBS，并参考已有成功脚本；
2. 询问用户 VASP/Wannier90 命令、队列、申请核数、实际 MPI ranks、内存、时限和环境设置；
3. 根据 INCAR 自动判断使用标准 VASP 还是 SOC/非共线版本；
4. 展示最终资源和启动命令；
5. 仅在用户明确确认后提交任务。

若 band 或 projection 任务出现可信的内存不足证据，agent 会询问用户选择减少实际 MPI ranks、增加内存申请，或同时采用两者，不会自动重提任务。

## 目录冲突保护

对新建的 `band/`、`proj/` 和 `wann/`：

- 目标目录不存在：使用原目录名；
- 目标目录已存在：新结果写入 `<name>_agent/`；
- 原目录和 `<name>_agent/` 都存在：停止并询问用户。

后续命令应读取 `stage_summary.json` 中的实际 `workdir`，并通过 `--band-dir`、`--proj-dir` 或 `--wann-dir` 继续传递，避免覆盖已有结果。

## 贝叶斯优化

当物理初始能窗未通过 LESA 时，可以优化：

- `dis_win_min`
- `dis_win_max`
- `dis_froz_min`
- `dis_froz_max`

每次生成新 batch 前，agent 必须询问用户希望并行运行几个候选。推荐值是 10，但实际数量完全服从用户指定：

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" \
  --root . \
  propose --batch-size <user-choice>
```

命令会生成 scheduler-independent trial 目录和 `optimization/batch-XXX.json`。每个候选由调度任务调用：

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" \
  --root . \
  run-trial \
  --trial-id <trial-id> \
  --wannier-command '<launcher> <wannier90.x> {{ seedname }}'
```

所有候选结束后：

```bash
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" --root . collect
python "$AUTO_MLWFS_SKILL_DIR/scripts/optimize_windows.py" --root . status
```

`skills/auto-mlwfs/examples/agi/optimization_phase2.yaml` 是针对 AgI 第一阶段中 `trial-0006` 撞到搜索边界后建立的第二阶段配置，不是通用默认配置。其他材料不应直接套用。

默认 `bo.max_evaluations` 为 `null`，不设置评估次数上限。需要控制计算预算时，应在项目配置中显式设置上限。使用自定义 `--config` 时，`propose`、`run-trial`、`collect` 和 `status` 应始终使用同一份配置。

## 机器可读输出

主要准备脚本会在实际工作目录写入 `stage_summary.json`，记录：

- 当前阶段和状态；
- 实际工作目录及是否使用 `_agent` 后缀；
- `NBANDS`、`NUM_WANN`、VASP 模式和物理能窗等关键参数；
- 尚待用户确认的调度器输入；
- 推荐的下一步动作。

BO 的 `optimization/batch-XXX.json` 还会记录候选参数、trial ID、工作目录和约束计数。上层 agent 应优先读取这些 JSON，而不是解析控制台文本。

## 配置与详细规则

内置配置相对于技能安装目录定位。修改参数时，建议先将 YAML 复制到计算项目，再通过 `--config` 指定：

- LESA 和 BO 的相对 `--config` 路径基于 `--root` 解析；DOS 的相对 `--config` 路径基于当前终端目录解析。
- YAML 中的计算数据路径相对于 `--root`，DOS 则相对于 `--workdir`。
- BO 的 `paths.fitting_config: null` 使用内置 LESA 配置；显式配置路径相对于 `--root`，也可使用绝对路径。

各阶段的详细说明：

- 投影规则：[`docs/projection_rules.md`](skills/auto-mlwfs/docs/projection_rules.md)
- band 工作流：[`docs/band_workflow.md`](skills/auto-mlwfs/docs/band_workflow.md)
- 调度器交互：[`docs/scheduler_guidance.md`](skills/auto-mlwfs/docs/scheduler_guidance.md)
- LESA 判据：[`docs/lesa_fitting.md`](skills/auto-mlwfs/docs/lesa_fitting.md)
- 成功标准：[`docs/success_criteria.md`](skills/auto-mlwfs/docs/success_criteria.md)
- 故障处理：[`docs/troubleshooting.md`](skills/auto-mlwfs/docs/troubleshooting.md)

## 使用建议与反馈

这是一个由实际材料计算流程逐步整理出来的初版 Skill。不同集群、VASP/Wannier90 编译方式和材料体系仍可能暴露新的适配问题。使用时请保留输入、日志、`stage_summary.json`、LESA 输出和触发问题的具体命令，以便后续复现和优化。

问题与建议请提交至 [GitHub Issues](https://github.com/sqwang2/auto-mlwfs/issues)。分享最小复现材料时，请移除集群凭据和受许可限制的文件。

## 开发与验证

在已准备好依赖的环境中，从仓库根目录运行：

```bash
python -m unittest discover -s tests -v
npx skills add . --list
```

现有测试覆盖技能安装可移植性和合成能带数据上的分析行为，CI 使用 Python 3.11 和 3.12。测试通过不等同于在所有材料、外部程序版本或集群上完成科学验证。

## 许可证与引用

本项目采用 [MIT 许可证](LICENSE)。外部科研软件保留各自的许可要求。

使用本软件开展研究时，可参考 [`CITATION.cff`](CITATION.cff) 引用项目，并按实际使用情况引用 VASP、Wannier90 及相关方法的原始文献。

---

**简体中文** | [English](README.en.md)
