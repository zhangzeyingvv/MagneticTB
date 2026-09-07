# MagneticTB

[English](README.md) | 简体中文

MagneticTB 用于构造受对称性约束的紧束缚哈密顿量，支持磁空间群、非磁空间群和自旋空间群。给定所选 Wyckoff 位置的对称性与轨道信息后，它可以生成对称性允许的哈密顿量，并提供对称操作、能带结构及相关紧束缚计算工具。

本仓库包含两个可以独立使用的实现：

- [`Mathematica/`](Mathematica/)：Wolfram Language 程序包、磁对称性数据、英文和简体中文文档以及示例 Notebook。
- [`Python-Rust-Web/`](Python-Rust-Web/)：Rust 计算核心、Python API、Web 界面、运行时数据和用户文档。

两个实现运行时互不依赖。

## 安装

### Mathematica

当前 Mathematica Paclet 版本为 `2.0.10`，需要 Wolfram Language 12.1 或更高版本。

从 Release 资源中下载 `MagneticTB-2.0.10.paclet`，然后在 Mathematica 内核中使用绝对路径安装：

```wl
PacletInstall["/absolute/path/to/MagneticTB-2.0.10.paclet"]
```

安装完成后退出并重新启动 Mathematica，然后在新内核中检查并载入程序包：

```wl
PacletFind["MagneticTB"]
Needs["MagneticTB`"]
```

要打开已安装的文档，请在 Mathematica 中选择 **帮助 > Wolfram 文档**，搜索 `MagneticTB`，然后在搜索结果中打开 **MagneticTB** 指南，即可浏览函数页面、教程和示例。也可以直接搜索 `init`、`initfromrep` 或 `symham` 等函数名，打开相应的参考页面。

更新、卸载、文档和示例的详细说明参见 [Mathematica README](Mathematica/README.zh-CN.md)。

### Python + Rust

Python 接口需要 CPython 3.9 或更高版本。从 Release 资源中下载与你的 Python 版本、操作系统和 CPU 架构匹配的 `.whl` 文件。Wheel 已包含编译后的 Rust 核心，普通用户不需要另外安装 Rust 或 Cargo。

在 macOS、Linux 或 WSL 中使用 `venv`：

```sh
python3 -m venv .venv
.venv/bin/python -m pip install /absolute/path/to/downloaded-wheel.whl
.venv/bin/python -c "import magnetictb; print(magnetictb.__version__)"
.venv/bin/magnetictb-web
```

在 Windows PowerShell 中使用 `venv`：

```powershell
py -3 -m venv .venv
.venv\Scripts\python.exe -m pip install C:\absolute\path\to\downloaded-wheel.whl
.venv\Scripts\python.exe -c "import magnetictb; print(magnetictb.__version__)"
.venv\Scripts\magnetictb-web.exe
```

推荐使用 Conda。请在已启用 Conda 的终端中运行；Windows 可以使用 Anaconda Prompt 或 Miniforge Prompt：

```sh
conda create -n magnetictb --override-channels -c conda-forge python=3.12 pip
conda activate magnetictb
python -m pip install /absolute/path/to/downloaded-wheel.whl
python -c "import magnetictb; print(magnetictb.__version__)"
magnetictb-web
```

当前程序包版本会输出 `0.1.0`。

启动 `magnetictb-web` 后，在浏览器中打开 <http://127.0.0.1:8000/>。交互式 API 文档位于 <http://127.0.0.1:8000/api/docs>。服务默认只监听本机；按 `Ctrl+C` 停止。

建模流程、API、精确输入和示例参见 [Python/Rust/Web 中文 README](Python-Rust-Web/README.zh-CN.md) 和 [Web 帮助中心](Python-Rust-Web/docs/user-guide/web/index.html)。

## 功能

- 构造受对称性约束的紧束缚哈密顿量。
- 使用磁空间群、非磁空间群和自旋空间群的对称性数据。
- 获得对称操作的矩阵表示。
- 按键层生成实空间和动量空间哈密顿量。
- 处理和分析能带结构及相关模型性质。
- 使用 Mathematica 接口，或使用由 Rust 计算核心支持的 Python API。

## 示例与文档

- Mathematica 示例：[`Mathematica/Examples/`](Mathematica/Examples/)
- Mathematica 双语帮助：[`Mathematica/Documentation/`](Mathematica/Documentation/)
- Python/Rust Web 帮助：[`Python-Rust-Web/docs/user-guide/web/`](Python-Rust-Web/docs/user-guide/web/index.html)

## 版本说明

版本按时间从新到旧排列。

### Python/Rust 0.1.0（2026-08-18）

- 加入 Rust 计算核心、Python API 和本地 Web 界面。

### Mathematica 2.0.10（2026-08-17）

- 扩充英文和简体中文文档、教程及快速入门材料。
- 加入实空间哈密顿量、薄片、表面格林函数、Berry 几何、Wilson loop、晶体、Brillouin 区和 k 路径工具，并改进能带可视化。

### Mathematica 2.0.0（2026-03-15）

- 利用线性代数和群表示论重写核心算法。
- 加入诱导表示模式，现在可以构造最小紧束缚模型。
- 加入 `initfromrep`，无需输入基函数，只需输入 site symmetry group 的表示即可构造紧束缚模型。
- 正式支持自旋空间群（SSG），包括 collinear、coplanar 和 non-coplanar 情形。
- 为 `symham` 加入分圆域精确零空间 `KernelMethod`。

### Mathematica 1.06（2025-12-11）

- 加入自旋空间群的测试支持。
- 加入磁空间群 Wyckoff 位置显示，以及磁性层群和杆群的对称操作数据。

### Mathematica 1.05（2024-12-04）

- 修复 `hop` 中的一个错误。
- 加入 `readHR`，用于导入 `wannier90_hr.dat` 文件。

### Mathematica 1.04（2024-06-21）

- 为 `symham` 加入 `CartesianCoordinates` 选项。
- 加入 `banddata`，用于生成 `band.dat` 文件。

### Mathematica 1.03（2023-02-17）

- 加入自动寻找空间群生成元的贪心算法，显著提高计算效率。

### Mathematica 1.02b（2023-02-14）

- 加入手工获得紧束缚参数的示例。
- 加入英文手册。

### Mathematica 1.02（2022-12-01）

- 使用可选的 `SpaceGroupIrep` 和 `MSGCorep` 程序包加入紧束缚能带余表示计算。
- 加入 `getMSGElemFromMSGCorep` 和 `getTBBandCorep`。

### Mathematica 1.01（2022-07-22）

- 修复自动酉化偶尔导致基函数顺序改变的问题。

### Mathematica 1.00c

- 加入中文手册。
- 修复单斜晶格矢量显示。

### Mathematica 1.00b

- 加入双磁空间群中电荷为 4 的 Weyl 点示例。

## 引用

如果 MagneticTB 对你的研究有帮助，请引用：

Z. Zhang, Z.-M. Yu, G.-B. Liu, and Y. Yao, “MagneticTB: A package for tight-binding model of magnetic and nonmagnetic materials,” *Computer Physics Communications* **270**, 108153 (2022).

- [期刊论文](https://www.sciencedirect.com/science/article/abs/pii/S0010465521002654)
- [arXiv:2105.09504](https://arxiv.org/abs/2105.09504)

## 许可证

MagneticTB 采用 GNU 通用公共许可证第 3 版（仅限该版本，
`GPL-3.0-only`）。详情见 [LICENSE](LICENSE)。

Copyright (C) 2021-2026 Zhang Zeying。
